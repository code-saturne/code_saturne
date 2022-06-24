/*============================================================================
 * Functions dedicated to the linear algebra settings and operations in case
 * of CDO face-based schemes with a monolithic velocity-pressure coupling
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <string.h>

#if defined(HAVE_OPENMP)
#include <omp.h>
#endif

#if defined(HAVE_PETSC)
#include <petscversion.h>
#include <petscksp.h>
#endif

/*----------------------------------------------------------------------------
 *  BFT headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_blas.h"
#include "cs_cdo_blas.h"
#include "cs_cdo_solve.h"
#include "cs_equation.h"
#include "cs_fp_exception.h"
#include "cs_matrix_default.h"
#include "cs_navsto_sles.h"
#include "cs_parall.h"
#include "cs_saddle_itsol.h"
#include "cs_timer.h"

#if defined(DEBUG) && !defined(NDEBUG)
#include "cs_dbg.h"
#endif

#if defined(HAVE_MUMPS)
#include "cs_sles_mumps.h"
#endif

#if defined(HAVE_PETSC)
#include "cs_sles_petsc.h"
#endif

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdofb_monolithic_sles.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
 * \file cs_cdofb_monolithic_sles.c
 *
 * \brief Functions dedicated to to the linear algebra settings and operations
 *        in case of CDO face-based schemes with a monolithic velocity-pressure
 *        coupling
 */

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_CDOFB_MONOLITHIC_SLES_DBG      0

/* GKB advanced settings */

#define CS_GKB_TRUNCATION_THRESHOLD       5

/* Block size for superblock algorithm */

#define CS_SBLOCK_BLOCK_SIZE 60

/* Cache line multiple, in cs_real_t units */

#define CS_CL  (CS_CL_SIZE/8)

/* Redefined the name of functions from cs_param_sles to get shorter names */

#define _petsc_cmd  cs_param_sles_petsc_cmd

/* This structure follow notations given in the article entitled
 * "An iterative generalized Golub-Kahan algorithm for problems in structural
 *  mechanics" by M. Arioli, C. Kruse, U. Ruede and N. Tardieu
 *
 * M space is isomorphic to the velocity space (size = 3.n_faces)
 * N space is isomorphic to the pressure space (size = n_cells)
 */

typedef struct {

  /* Value of the grad-div coefficient */

  cs_real_t                gamma;

  /* Size of spaces */

  cs_lnum_t                n_u_dofs; /* Size of the space M */
  cs_lnum_t                n_p_dofs; /* Size of the space N */

  /* Vector transformation */

  cs_real_t               *b_tilda;  /* Modified RHS */
  cs_real_t               *u_tilda;  /* Modified velocity unknown */

  /* Auxiliary vectors */

  cs_real_t               *q;        /* vector iterates in space N */
  cs_real_t               *d;        /* vector iterates in space N */
  cs_real_t               *d__v;     /* buffer in space N */
  cs_real_t               *dt_q;     /* buffer in space M */
  cs_real_t               *m__v;     /* vector iterates in space M */
  cs_real_t               *v;        /* vector iterates in space M */

  /* Orthogonalization coefficients */

  cs_real_t                alpha;
  cs_real_t                beta;
  cs_real_t                zeta;

  /* Store z_size zeta coefficients */

  int                      z_size;
  cs_real_t               *zeta_array;
  cs_real_t                zeta_square_sum;

  cs_iter_algo_t          *algo;     /* Information related to the convergence
                                        of the algorithm */

} cs_gkb_builder_t;

/* This structure is used to manage the Uzawa algorithm and its variants
 *
 * U space is isomorphic to the velocity space (size = 3.n_faces)
 * P space is isomorphic to the pressure space (size = n_cells)
 */

typedef struct {

  /* Value of the grad-div coefficient */

  cs_real_t               gamma;

  /* Size of spaces */

  cs_lnum_t               n_u_dofs; /* Size of the space U */
  cs_lnum_t               n_p_dofs; /* Size of the space P */

  /* Vector transformation */

  cs_real_t              *b_tilda;  /* Modified RHS (size U) */

  /* Auxiliary scaling coefficient */

  cs_real_t               alpha;

  /* Auxiliary vectors */

  cs_real_t              *inv_mp;   /* reciprocal of the pressure mass matrix */
  cs_real_t              *res_p;    /* buffer in space P */
  cs_real_t              *d__v;     /* buffer in space P */
  cs_real_t              *gk;       /* buffer in space P */
  cs_real_t              *dzk;      /* buffer in space U */
  cs_real_t              *rhs;      /* buffer in space U */

  cs_iter_algo_t         *algo;     /* Information related to the convergence
                                       of the algorithm */

} cs_uza_builder_t;

/* Context structure for more complex PETSc configurations */

typedef struct {

  const cs_navsto_param_t      *nsp;
  const cs_cdofb_monolithic_t  *sc;

} cs_cdofb_monolithic_petsc_context_t;

/*============================================================================
 * Private variables
 *============================================================================*/

/* Pointer to shared structures */

static const cs_cdo_connect_t       *cs_shared_connect;
static const cs_cdo_quantities_t    *cs_shared_quant;

static cs_cdofb_monolithic_petsc_context_t   *_petsc_hook_context = NULL;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute a norm for a scalar-valued cell-based array "a"
 *        The parallel synchronization is performed inside this function
 *
 * \param[in]    a    array of size n_cells
 *
 * \return the computed norm
 */
/*----------------------------------------------------------------------------*/

static inline double
_get_cbscal_norm(cs_real_t  *a)
{
  double norm2 = cs_dot_xx(cs_shared_quant->n_cells, a);

  cs_parall_sum(1, CS_DOUBLE, &norm2);
  assert(norm2 > -DBL_MIN);

  return sqrt(norm2);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute a norm for a face-based "a" v with 3*n_faces elements
 *        The parallel synchronization is performed inside this function
 *
 * \param[in]    a    array of size 3*n_faces
 *
 * \return the computed norm
 */
/*----------------------------------------------------------------------------*/

static inline double
_get_fbvect_norm(cs_real_t  *a)
{
  double norm2 = cs_cdo_blas_square_norm_pfvp(a);

  assert(norm2 > -DBL_MIN);
  return sqrt(norm2);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute array index bounds for a local thread.
 *
 * When called inside an OpenMP parallel section, this will return the
 * start an past-the-end indexes for the array range assigned to that thread.
 * In other cases, the start index is 1, and the past-the-end index is n;
 *
 * \param[in]   n     size of array
 * \param[out]  s_id  start index for the current thread
 * \param[out]  e_id  past-the-end index for the current thread
 */
/*----------------------------------------------------------------------------*/

static inline void
_thread_range(cs_lnum_t   n,
              cs_lnum_t  *s_id,
              cs_lnum_t  *e_id)
{
#if defined(HAVE_OPENMP)
  int t_id = omp_get_thread_num();
  int n_t = omp_get_num_threads();
  cs_lnum_t t_n = (n + n_t - 1) / n_t;
  *s_id =  t_id    * t_n;
  *e_id = (t_id+1) * t_n;
  *s_id = cs_align(*s_id, CS_CL);
  *e_id = cs_align(*e_id, CS_CL);
  if (*e_id > n) *e_id = n;
#else
  *s_id = 0;
  *e_id = n;
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Dot product between two arrays on face unknowns.
 *         One assumes that input arrays are in a "scattered" distribution
 *         So the size should be 3*n_faces.
 *
 * \param[in]       rset   pointer to a range_set structure (synchro. op.)
 * \param[in]       size   size of arrays
 * \param[in, out]  x      first array
 * \param[in, out]  y      second array
 *
 * \return the computed value
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
_face_gdot(const cs_range_set_t   *rset,
           cs_lnum_t               size,
           cs_real_t               x[],
           cs_real_t               y[])
{
  CS_UNUSED(size); /* Avoid a compilation warning in during compilation */
  assert(size == rset->n_elts[1]);
  assert(size == 3*cs_shared_quant->n_faces);

  /* x and y are scattered arrays. One assumes that values are synchronized
     across ranks (for instance by using a cs_interface_set_sum()) */

  cs_range_set_gather(rset,
                      CS_REAL_TYPE,/* type */
                      1,           /* stride (treated as scalar up to now) */
                      x,           /* in: size = n_sles_scatter_elts */
                      x);          /* out: size = n_sles_gather_elts */

  cs_range_set_gather(rset,
                      CS_REAL_TYPE,/* type */
                      1,           /* stride (treated as scalar up to now) */
                      y,           /* in: size = n_sles_scatter_elts */
                      y);          /* out: size = n_sles_gather_elts */

  cs_real_t  result = cs_gdot(rset->n_elts[0], x, y);

  cs_range_set_scatter(rset,
                       CS_REAL_TYPE,
                       1,
                       x,
                       x);
  cs_range_set_scatter(rset,
                       CS_REAL_TYPE,
                       1,
                       y,
                       y);

  return result;
}

#if defined(HAVE_PETSC)
/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize the cs_cdofb_monolithic_petsc_context_t
 *
 * \param[in]  nsp      pointer to the set of parameters related to NavSto
 * \param[in]  sc       scheme context related to CDO-Fb monolithic schemes
 */
/*----------------------------------------------------------------------------*/

static inline void
_initialize_petsc_hook_context(cs_navsto_param_t      *nsp,
                               cs_cdofb_monolithic_t  *sc)
{
  /* Initialization must be called before setting options;
     it does not need to be called before calling
     cs_sles_petsc_define(), as this is handled automatically. */

  if (_petsc_hook_context == NULL)
    BFT_MALLOC(_petsc_hook_context, 1, cs_cdofb_monolithic_petsc_context_t);

  _petsc_hook_context->nsp = nsp;
  _petsc_hook_context->sc = sc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the default options for a PCGAMG type in PETSc
 *
 * \param[in]  prefix        optional prefix
 * \param[in]  system_size   size of the linear system
 * \param[in]  amg_type      type of AMG preconditioner
 * \param[in]  is_sym        system to solve is symmetric ?
 * \param[in]  smooth_lvl    level of smoothing (0: light)
 */
/*----------------------------------------------------------------------------*/

static void
_set_gamg_pc(const char            prefix[],
             cs_gnum_t             system_size,
             cs_param_amg_type_t   amg_type,
             bool                  is_sym,
             int                   smooth_lvl)
{
  /* Estimate the number of levels */

  double  _n_levels = ceil(0.6*(log(system_size) - 5));
  int n_levels = 1, max_levels = 14;
  if (_n_levels > 1)
    n_levels = (int)_n_levels;

  /* Need to add a prefix */

  bool  use_pre = (prefix != NULL) ? true : false;
  if (use_pre)
    if (strlen(prefix) < 1) use_pre = false;

  /* Set the type of cycles (V or W) */

  switch(amg_type) {

  case CS_PARAM_AMG_HYPRE_BOOMER_V:
  case CS_PARAM_AMG_PETSC_GAMG_V:
    _petsc_cmd(use_pre, prefix, "pc_mg_cycle_type", "v");
    n_levels = CS_MIN(n_levels, max_levels) + 1;
    break;
  case CS_PARAM_AMG_HYPRE_BOOMER_W:
  case CS_PARAM_AMG_PETSC_GAMG_W:
    _petsc_cmd(use_pre, prefix, "pc_mg_cycle_type", "w");
    n_levels = CS_MIN(n_levels, max_levels);
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Incompatible type of AMG for PETSc.\n",
              __func__);
    break;
  }

  char  string_n_levels[4];
  sprintf(string_n_levels, "%2d", n_levels);
  _petsc_cmd(use_pre, prefix, "pc_mg_levels", string_n_levels);

  /* Symmetrize the graph before computing the aggregation. Some algorithms
   * require the graph be symmetric (default=false) */

  _petsc_cmd(use_pre, prefix, "pc_gamg_sym_graph", "true");

  /* Remark: -pc_gamg_reuse_interpolation
   *
   * Reuse prolongation when rebuilding algebraic multigrid
   * preconditioner. This may negatively affect the convergence rate of the
   * method on new matrices if the matrix entries change a great deal, but
   * allows rebuilding the preconditioner quicker. (default=false)
   */

  /* Remark: -pc_gamg_square_graph
   *
   * Squaring the graph increases the rate of coarsening (aggressive
   * coarsening) and thereby reduces the complexity of the coarse grids, and
   * generally results in slower solver converge rates. Reducing coarse grid
   * complexity reduced the complexity of Galerkin coarse grid construction
   * considerably. (default = 1)
   *
   * Remark: -pc_gamg_threshold
   *
   * Increasing the threshold decreases the rate of coarsening. Conversely
   * reducing the threshold increases the rate of coarsening (aggressive
   * coarsening) and thereby reduces the complexity of the coarse grids, and
   * generally results in slower solver converge rates. Reducing coarse grid
   * complexity reduced the complexity of Galerkin coarse grid construction
   * considerably. Before coarsening or aggregating the graph, GAMG removes
   * small values from the graph with this threshold, and thus reducing the
   * coupling in the graph and a different (perhaps better) coarser set of
   * points. (default=0.0) */

  _petsc_cmd(use_pre, prefix, "mg_levels_ksp_norm_type", "none");
  _petsc_cmd(use_pre, prefix, "pc_gamg_coarse_eq_limit", "100");

  /* In parallel computing, migrate data to another rank if the grid has less
     than 200 rows */

  if (cs_glob_n_ranks > 1) {
    _petsc_cmd(use_pre, prefix, "pc_gamg_repartition", "true");
    _petsc_cmd(use_pre, prefix, "pc_gamg_process_eq_limit", "200");
  }

  /* More efficient sparse direct solver */

  if (cs_glob_n_ranks == 1) {
    _petsc_cmd(use_pre, prefix, "mg_coarse_ksp_type", "preonly");
    _petsc_cmd(use_pre, prefix, "mg_coarse_pc_type", "tfs");
  }

  if (is_sym) {

    /* Number of smoothing steps to use with smooth aggregation (default=1) */

    _petsc_cmd(use_pre, prefix, "pc_gamg_agg_nsmooths", "1");
    _petsc_cmd(use_pre, prefix, "pc_gamg_reuse_interpolation", "false");
    _petsc_cmd(use_pre, prefix, "pc_gamg_esteig_ksp_type", "cg");

    /* PCMG settings (options shared with PCGAMG) */

    _petsc_cmd(use_pre, prefix, "pc_gamg_threshold", "0.10");
    _petsc_cmd(use_pre, prefix, "pc_gamg_square_graph", "2");

    /* Apply one Richardson relaxation (scaling = 1.0 -- the default value) */

    _petsc_cmd(use_pre, prefix, "mg_levels_ksp_max_it", "1");
    _petsc_cmd(use_pre, prefix, "mg_levels_ksp_type", "richardson");
    _petsc_cmd(use_pre, prefix, "mg_levels_ksp_richardson_scale", "1.0");

    /* Set the up/down smoothers */

    if (smooth_lvl == 0) {

      _petsc_cmd(use_pre, prefix, "mg_levels_pc_type", "sor");

    }
    else if (smooth_lvl == 1) {

      _petsc_cmd(use_pre, prefix, "mg_levels_pc_type", "sor");
      _petsc_cmd(use_pre, prefix, "mg_levels_pc_sor_lits", "2");

    }
    else {

      /* Each Richardson step is completed with a ILU(0)-GMRES call */

      _petsc_cmd(use_pre, prefix, "mg_levels_pc_type", "bjacobi");
      _petsc_cmd(use_pre, prefix, "mg_levels_pc_bjacobi_blocks", "1");
      _petsc_cmd(use_pre, prefix, "mg_levels_sub_ksp_type", "gmres");
      _petsc_cmd(use_pre, prefix, "mg_levels_sub_ksp_gmres_restart", "5");
      _petsc_cmd(use_pre, prefix, "mg_levels_sub_ksp_max_it", "5");
      _petsc_cmd(use_pre, prefix, "mg_levels_sub_ksp_rtol", "1e-5");
      _petsc_cmd(use_pre, prefix, "mg_levels_sub_pc_type", "sor");

    }

  }
  else { /* Not a symmetric system */

    /* Number of smoothing steps to use with smooth aggregation (default=1) */

    _petsc_cmd(use_pre, prefix, "pc_gamg_agg_nsmooths", "0");
    _petsc_cmd(use_pre, prefix, "pc_gamg_reuse_interpolation", "false");
    _petsc_cmd(use_pre, prefix, "pc_gamg_square_graph", "0");
    _petsc_cmd(use_pre, prefix, "pc_gamg_threshold", "0.06");

    /* Set the up/down smoothers
     * Apply one Richardson relaxation (scaling = 1.0 -- the default value) */

    _petsc_cmd(use_pre, prefix, "mg_levels_ksp_max_it", "1");
    _petsc_cmd(use_pre, prefix, "mg_levels_ksp_type", "richardson");
    _petsc_cmd(use_pre, prefix, "mg_levels_ksp_richardson_scale", "1.0");

    if (smooth_lvl == 0) {

      _petsc_cmd(use_pre, prefix, "mg_levels_pc_type", "sor");

    }
    else if (smooth_lvl == 1) {

      _petsc_cmd(use_pre, prefix, "mg_levels_pc_type", "bjacobi");
      _petsc_cmd(use_pre, prefix, "mg_levels_pc_bjacobi_blocks", "1");
      _petsc_cmd(use_pre, prefix, "mg_levels_sub_ksp_type", "preonly");
      _petsc_cmd(use_pre, prefix, "mg_levels_sub_pc_type", "ilu");
      _petsc_cmd(use_pre, prefix, "mg_levels_sub_pc_factor_levels", "0");

    }
    else {

      /* Each Richardson step is completed with a ILU(0)-GMRES call */

      _petsc_cmd(use_pre, prefix, "mg_levels_pc_type", "bjacobi");
      _petsc_cmd(use_pre, prefix, "mg_levels_pc_bjacobi_blocks", "1");
      _petsc_cmd(use_pre, prefix, "mg_levels_sub_ksp_type", "gmres");
      _petsc_cmd(use_pre, prefix, "mg_levels_sub_ksp_gmres_restart", "10");
      _petsc_cmd(use_pre, prefix, "mg_levels_sub_ksp_max_it", "10");
      _petsc_cmd(use_pre, prefix, "mg_levels_sub_ksp_rtol", "1e-5");
      _petsc_cmd(use_pre, prefix, "mg_levels_sub_ksp_norm_type", "none");
      _petsc_cmd(use_pre, prefix, "mg_levels_sub_pc_type", "ilu");
      _petsc_cmd(use_pre, prefix, "mg_levels_sub_pc_factor_levels", "0");

    } /* Light smoothing ? */

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the way to normalize the residual vector
 *
 * \param[in]       norm_type   type of normalization
 * \param[in, out]  ksp         pointer to a PETSc KSP structure
 */
/*----------------------------------------------------------------------------*/

static inline void
_set_residual_normalization(cs_param_resnorm_type_t    norm_type,
                            KSP                        ksp)
{
  switch (norm_type) {

  case CS_PARAM_RESNORM_NORM2_RHS: /* Try to have "true" norm */
    KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED);
    break;
  case CS_PARAM_RESNORM_NONE:
    KSPSetNormType(ksp, KSP_NORM_NONE);
    break;
  default:
    KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED);
    break;

  }
}
#endif  /* HAVE_PETSC */

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the matrix for an approximation of the Schur complement based
 *         on the inverse of the diagonal of the velocity block
 *
 * \param[in]      nsp          pointer to a cs_navsto_param_t structure
 * \param[in]      a            (MSR) matrix for the velocity block
 * \param[in]      rset         pointer to the associated range set structure
 * \param[in, out] uza          structure to manage the Uzawa algorithm
 * \param[out]     p_diag_smat  diagonal coefficients for the Schur matrix
 * \param[out]     p_xtra_smat  extra-diagonal coefficients for the Schur matrix
 *
 * \return a pointer to the computed matrix
 */
/*----------------------------------------------------------------------------*/

static cs_matrix_t *
_diag_schur_approximation(const cs_navsto_param_t   *nsp,
                          const cs_matrix_t         *a,
                          const cs_range_set_t      *rset,
                          cs_uza_builder_t          *uza,
                          cs_real_t                **p_diag_smat,
                          cs_real_t                **p_xtra_smat)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_lnum_t  n_cells = m->n_cells;
  const cs_lnum_t  n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t  n_i_faces = m->n_i_faces;
  const cs_lnum_t  n_b_faces = m->n_b_faces;
  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  /* Synchronize the diagonal values for A */

  const cs_real_t  *diagA = NULL;
  cs_real_t  *_diagA = NULL;

  if (cs_glob_n_ranks > 1) {

    size_t  size = 3*quant->n_faces;
    BFT_MALLOC(_diagA, size, cs_real_t);
    cs_range_set_scatter(rset,
                         CS_REAL_TYPE,
                         1,     /* treated as scalar-valued up to now */
                         cs_matrix_get_diagonal(a), /* gathered view */
                         _diagA);

    diagA = _diagA; /* scatter view (synchronized)*/

  }
  else {

    diagA = cs_matrix_get_diagonal(a);
    assert(m->periodicity == NULL); /* TODO */

  }

  /* Native format for the Schur approximation matrix */

  cs_real_t   *diag_smat = NULL;
  cs_real_t   *xtra_smat = NULL;

  BFT_MALLOC(diag_smat, n_cells_ext, cs_real_t);
  BFT_MALLOC(xtra_smat, 2*n_i_faces, cs_real_t);

  memset(diag_smat, 0, n_cells_ext*sizeof(cs_real_t));
  memset(xtra_smat, 0, 2*n_i_faces*sizeof(cs_real_t));

  /* Add diagonal and extra-diagonal contributions from interior faces */

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {

    const cs_real_t  *a_ff = diagA + 3*f_id;
    const cs_nvec3_t  nvf = cs_quant_set_face_nvec(f_id, quant);

    double  contrib = 0;
    for (int k = 0; k < 3; k++)
      contrib += 1/a_ff[k]*nvf.unitv[k]*nvf.unitv[k];
    contrib *= -nvf.meas*nvf.meas;

    /* Extra-diagonal contribution. This is scanned by the i_face_cells mesh
       adjacency */

    cs_real_t  *_xtra_smat = xtra_smat + 2*f_id;
    _xtra_smat[0] = _xtra_smat[1] = contrib;

    /* Diagonal contributions */

    cs_lnum_t cell_i = i_face_cells[f_id][0];
    cs_lnum_t cell_j = i_face_cells[f_id][1];

    diag_smat[cell_i] -= contrib;
    diag_smat[cell_j] -= contrib;

  } /* Loop on interior faces */

  /* Add diagonal contributions from border faces*/

  const cs_real_t  *diagA_shift = diagA + 3*n_i_faces;
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

    const cs_real_t  *a_ff = diagA_shift + 3*f_id;

    cs_nvec3_t  nvf;
    cs_nvec3(quant->b_face_normal + 3*f_id, &nvf);

    double  contrib = 0;
    for (int k = 0; k < 3; k++)
      contrib += 1/a_ff[k]*nvf.unitv[k]*nvf.unitv[k];
    contrib *= nvf.meas*nvf.meas;

    /* Diagonal contributions */

    diag_smat[b_face_cells[f_id]] += contrib;

  } /* Loop on border faces */

  /* Return the associated matrix */

  /* One assumes a non-symmetric matrix even if in most (all?) cases the matrix
     should be symmetric */

  cs_matrix_t  *smat = NULL;

  if (nsp->sles_param->schur_sles_param->solver_class ==
      CS_PARAM_SLES_CLASS_HYPRE)
    smat = cs_matrix_external("HYPRE_ParCSR",
                              false, /* symmetry */
                              1, 1);
  else
    smat = cs_matrix_msr(false, /* symmetry */
                         1, 1);

  cs_matrix_set_coefficients(smat, false, /* symmetry */
                             1, 1,
                             n_i_faces, i_face_cells,
                             diag_smat, xtra_smat);

  /* Return arrays (to be freed when the algorithm is converged) */

  *p_diag_smat = diag_smat;
  *p_xtra_smat = xtra_smat;

  BFT_FREE(_diagA);

  return smat;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the matrix for an approximation of the Schur complement based
 *         on the inverse of the sum of the absolute values of the velocity
 *         block
 *
 * \param[in]      nsp          pointer to a cs_navsto_param_t structure
 * \param[in]      eqp          pointer to the set of equation parameters
 * \param[in]      a            (MSR) matrix for the velocity block
 * \param[in]      rset         pointer to a range set structure
 * \param[in, out] slesp        pointer to a set of parameters to drive the SLES
 * \param[in, out] msles        structure to manage the monolithic SLES
 * \param[in, out] uza          structure to manage the Uzawa algorithm
 * \param[out]     p_diag_smat  diagonal coefficients for the Schur matrix
 * \param[out]     p_xtra_smat  extra-diagonal coefficients for the Schur matrix
 *
 * \return a pointer to the computed matrix
 */
/*----------------------------------------------------------------------------*/

static cs_matrix_t *
_invlumped_schur_approximation(const cs_navsto_param_t     *nsp,
                               const cs_equation_param_t   *eqp,
                               const cs_matrix_t           *a,
                               const cs_range_set_t        *rset,
                               cs_param_sles_t             *slesp,
                               cs_cdofb_monolithic_sles_t  *msles,
                               cs_uza_builder_t            *uza,
                               cs_real_t                  **p_diag_smat,
                               cs_real_t                  **p_xtra_smat)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_lnum_t  n_cells = m->n_cells;
  const cs_lnum_t  n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t  n_i_faces = m->n_i_faces;
  const cs_lnum_t  n_b_faces = m->n_b_faces;
  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  /* Compute A^-1 lumped. Consider a rhs with only one values. */

  for (cs_lnum_t i = 0; i < uza->n_u_dofs; i++)
    uza->rhs[i] = 1;

  cs_real_t  *invA_lumped = NULL;
  BFT_MALLOC(invA_lumped, uza->n_u_dofs, cs_real_t);
  memset(invA_lumped, 0, sizeof(cs_real_t)*uza->n_u_dofs);

  /* Modify the tolerance. Only a coarse approximation is needed */

  char  *init_system_name = slesp->name;
  double  init_eps = slesp->eps;
  int  init_max_iter = slesp->n_max_iter;

  char  *system_name = NULL;
  BFT_MALLOC(system_name, strlen(eqp->name) + strlen(":inv_lumped") + 1, char);
  sprintf(system_name, "%s:inv_lumped", eqp->name);

  slesp->name = system_name;
  slesp->eps = 1e-2;  /* Only a coarse approximation is needed */
  slesp->n_max_iter = 10;

  cs_param_sles_update_cvg_settings(true, slesp); /* use the field id */

  uza->algo->n_inner_iter
    += (uza->algo->last_inner_iter =
        cs_cdo_solve_scalar_system(uza->n_u_dofs,
                                   slesp,
                                   a,
                                   rset,
                                   1,     /* no normalization */
                                   false, /* rhs_redux --> already done */
                                   msles->sles,
                                   invA_lumped,
                                   uza->rhs));

  /* Set back the initial parameters */

  slesp->name = init_system_name;
  slesp->eps = init_eps;
  slesp->n_max_iter = init_max_iter;

  cs_param_sles_update_cvg_settings(true, slesp); /* use the field id */

  /* Partial memory free */

  BFT_FREE(system_name);

  /* Native format for the Schur approximation matrix */

  cs_real_t   *diag_smat = NULL;
  cs_real_t   *xtra_smat = NULL;

  BFT_MALLOC(diag_smat, n_cells_ext, cs_real_t);
  BFT_MALLOC(xtra_smat, 2*n_i_faces, cs_real_t);

  memset(diag_smat, 0, n_cells_ext*sizeof(cs_real_t));
  memset(xtra_smat, 0, 2*n_i_faces*sizeof(cs_real_t));

  /* Add diagonal and extra-diagonal contributions from interior faces */

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {

    const cs_real_t  *a_ff = invA_lumped + 3*f_id;
    const cs_nvec3_t  nvf = cs_quant_set_face_nvec(f_id, quant);

    double  contrib = 0;
    for (int k = 0; k < 3; k++)
      contrib += a_ff[k]*nvf.unitv[k]*nvf.unitv[k];
    contrib *= -nvf.meas*nvf.meas;

    /* Extra-diagonal contribution. This is scanned by the i_face_cells mesh
       adjacency */

    cs_real_t  *_xtra_smat = xtra_smat + 2*f_id;
    _xtra_smat[0] = _xtra_smat[1] = contrib;

    /* Diagonal contributions */

    cs_lnum_t cell_i = i_face_cells[f_id][0];
    cs_lnum_t cell_j = i_face_cells[f_id][1];

    diag_smat[cell_i] -= contrib;
    diag_smat[cell_j] -= contrib;

  } /* Loop on interior faces */

  /* Add diagonal contributions from border faces */

  cs_real_t  *diagA_shift = invA_lumped + 3*n_i_faces;
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

    const cs_real_t  *a_ff = diagA_shift + 3*f_id;

    cs_nvec3_t  nvf;
    cs_nvec3(quant->b_face_normal + 3*f_id, &nvf);

    double  contrib = 0;
    for (int k = 0; k < 3; k++)
      contrib += a_ff[k]*nvf.unitv[k]*nvf.unitv[k];
    contrib *= nvf.meas*nvf.meas;

    /* Diagonal contributions */

    diag_smat[b_face_cells[f_id]] += contrib;

  } /* Loop on border faces */

  /* Return the associated matrix */

  /* One assumes a non-symmetric matrix even if in most (all?) cases the matrix
     should be symmetric */

  cs_matrix_t  *smat = NULL;

  if (nsp->sles_param->schur_sles_param->solver_class ==
      CS_PARAM_SLES_CLASS_HYPRE)
    smat = cs_matrix_external("HYPRE_ParCSR",
                              false, /* symmetry */
                              1, 1);
  else
    smat = cs_matrix_msr(false, /* symmetry */
                         1, 1);

  cs_matrix_set_coefficients(smat, false, /* symmetry */
                             1, 1,
                             n_i_faces, i_face_cells,
                             diag_smat, xtra_smat);

  /* Return arrays (to be freed when the algorithm is converged) */

  *p_diag_smat = diag_smat;
  *p_xtra_smat = xtra_smat;

  BFT_FREE(invA_lumped);

  return smat;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the matrix for an approximation of the Schur complement based
 *        on the inverse of the sum of the absolute values of the velocity
 *        block
 *
 * \param[in]      nsp     pointer to a cs_navsto_param_t structure
 * \param[in]      ssys    pointer to a saddle-point system structure
 * \param[in, out] sbp     pointer to a cs_saddle_block_precond_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_diag_schur_sbp(const cs_navsto_param_t       *nsp,
                const cs_saddle_system_t      *ssys,
                cs_saddle_block_precond_t     *sbp)
{
  CS_UNUSED(nsp);

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_lnum_t  n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t  n_i_faces = m->n_i_faces;
  const cs_lnum_t  n_b_faces = m->n_b_faces;
  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  const cs_lnum_t  b11_size = ssys->x1_size;

  /* Synchronize the diagonal values for the block m11 */

  const cs_matrix_t  *m11 = ssys->m11_matrices[0];
  const cs_lnum_t  n_rows = cs_matrix_get_n_rows(m11);
  const cs_real_t  *diag_m11 = cs_matrix_get_diagonal(m11);

  cs_real_t  *inv_diag = NULL;
  BFT_MALLOC(inv_diag, CS_MAX(b11_size, n_rows), cs_real_t);

  /*  Operation in gather view (the default view for a matrix) */

  for (cs_lnum_t i1 = 0; i1 < n_rows; i1++)
    inv_diag[i1] = 1./diag_m11[i1];

  cs_range_set_scatter(ssys->rset,
                       CS_REAL_TYPE, 1, /* treated as scalar-valued up to now */
                       inv_diag,        /* gathered view */
                       inv_diag);       /* scatter view */

  /* Native format for the Schur approximation matrix */

  cs_real_t   *diag_smat = NULL;
  cs_real_t   *xtra_smat = NULL;

  BFT_MALLOC(diag_smat, n_cells_ext, cs_real_t);
  BFT_MALLOC(xtra_smat, 2*n_i_faces, cs_real_t);

  memset(diag_smat, 0, n_cells_ext*sizeof(cs_real_t));
  memset(xtra_smat, 0, 2*n_i_faces*sizeof(cs_real_t));

  /* Add diagonal and extra-diagonal contributions from interior faces */

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {

    const cs_real_t  *ia_ff = inv_diag + 3*f_id;
    const cs_nvec3_t  nvf = cs_quant_set_face_nvec(f_id, quant);

    double  contrib = 0;
    for (int k = 0; k < 3; k++)
      contrib += ia_ff[k]*nvf.unitv[k]*nvf.unitv[k];
    contrib *= -nvf.meas*nvf.meas;

    /* Extra-diagonal contribution. This is scanned by the i_face_cells mesh
       adjacency */

    cs_real_t  *_xtra_smat = xtra_smat + 2*f_id;
    _xtra_smat[0] = _xtra_smat[1] = contrib;

    /* Diagonal contributions */

    cs_lnum_t cell_i = i_face_cells[f_id][0];
    cs_lnum_t cell_j = i_face_cells[f_id][1];

    diag_smat[cell_i] -= contrib;
    diag_smat[cell_j] -= contrib;

  } /* Loop on interior faces */

  /* Add diagonal contributions from border faces*/

  const cs_real_t  *diag_shift = inv_diag + 3*n_i_faces;
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

    const cs_real_t  *ia_ff = diag_shift + 3*f_id;

    cs_nvec3_t  nvf;
    cs_nvec3(quant->b_face_normal + 3*f_id, &nvf);

    double  contrib = 0;
    for (int k = 0; k < 3; k++)
      contrib += ia_ff[k]*nvf.unitv[k]*nvf.unitv[k];
    contrib *= nvf.meas*nvf.meas;

    /* Diagonal contributions */

    diag_smat[b_face_cells[f_id]] += contrib;

  } /* Loop on border faces */

  /* Return the associated matrix */

  /* One assumes a non-symmetric matrix even if in most (all?) cases the matrix
     should be symmetric */

  if (nsp->sles_param->schur_sles_param->solver_class ==
      CS_PARAM_SLES_CLASS_HYPRE)
    sbp->schur_matrix = cs_matrix_external("HYPRE_ParCSR",
                                           false, /* symmetry */
                                           1, 1);
  else
    sbp->schur_matrix = cs_matrix_msr(false, /* symmetry */
                                      1, 1);

  cs_matrix_set_coefficients(sbp->schur_matrix, false, /* symmetry */
                             1, 1,
                             n_i_faces, i_face_cells,
                             diag_smat, xtra_smat);

  /* Return arrays (to be freed when the algorithm is converged) */

  sbp->schur_diag = diag_smat;
  sbp->schur_xtra = xtra_smat;
  sbp->m11_inv_diag = inv_diag;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a scaled mass matrix (on the pressure space) and a scaling
 *        coefficient for the compatible Laplacian
 *
 * \param[in]      nsp     pointer to a cs_navsto_param_t structure
 * \param[in]      ssys    pointer to a saddle-point system structure
 * \param[in, out] sbp     pointer to a cs_saddle_block_precond_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_scaled_mass_sbp(const cs_navsto_param_t       *nsp,
                 const cs_saddle_system_t      *ssys,
                 cs_saddle_block_precond_t     *sbp)
{
  CS_UNUSED(sbp);
  CS_UNUSED(ssys);

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_time_step_t  *ts = cs_glob_time_step;
  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_lnum_t  n_cells = m->n_cells;

  assert(ssys->x2_size == n_cells);
  BFT_MALLOC(sbp->mass22_diag, n_cells, cs_real_t);

  /* Compute scaling coefficients */

  if (nsp->turbulence->model->iturb == CS_TURB_NONE) {

    const cs_real_t  visc_val = nsp->lam_viscosity->ref_value;
#   pragma omp parallel for if (n_cells > CS_THR_MIN)
    for (cs_lnum_t i2 = 0; i2 < n_cells; i2++)
      sbp->mass22_diag[i2] = visc_val/quant->cell_vol[i2];

  }
  else {

    cs_property_eval_at_cells(ts->t_cur, nsp->tot_viscosity, sbp->mass22_diag);
#   pragma omp parallel for if (n_cells > CS_THR_MIN)
    for (cs_lnum_t i2 = 0; i2 < n_cells; i2++)
      sbp->mass22_diag[i2] /= quant->cell_vol[i2];

  }

  const cs_real_t  rho0 = nsp->mass_density->ref_value;
  cs_real_t  alpha = 1/ts->dt[0];
  if (nsp->model_flag & CS_NAVSTO_MODEL_STEADY)
    alpha = 0.01*nsp->lam_viscosity->ref_value;
  sbp->schur_scaling = rho0*alpha;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the matrix for an approximation of the Schur complement based
 *        on the inverse of the sum of the absolute values of the velocity
 *        block
 *
 * \param[in]      nsp     pointer to a cs_navsto_param_t structure
 * \param[in]      ssys    pointer to a saddle-point system structure
 * \param[in, out] sbp     pointer to a cs_saddle_block_precond_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_elman_schur_sbp(const cs_navsto_param_t       *nsp,
                 const cs_saddle_system_t      *ssys,
                 cs_saddle_block_precond_t     *sbp)
{
  CS_UNUSED(ssys);

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_lnum_t  n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t  n_i_faces = m->n_i_faces;
  const cs_lnum_t  n_b_faces = m->n_b_faces;
  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  /* Native format for the Schur approximation matrix */

  cs_real_t   *diag_smat = NULL;
  cs_real_t   *xtra_smat = NULL;

  BFT_MALLOC(diag_smat, n_cells_ext, cs_real_t);
  BFT_MALLOC(xtra_smat, 2*n_i_faces, cs_real_t);

  memset(diag_smat, 0, n_cells_ext*sizeof(cs_real_t));
  memset(xtra_smat, 0, 2*n_i_faces*sizeof(cs_real_t));

  /* Add diagonal and extra-diagonal contributions from interior faces */

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {

    const cs_nvec3_t  nvf = cs_quant_set_face_nvec(f_id, quant);

    double  contrib = 0;
    for (int k = 0; k < 3; k++)
      contrib += nvf.unitv[k]*nvf.unitv[k];
    contrib *= -nvf.meas*nvf.meas;

    /* Extra-diagonal contribution. This is scanned by the i_face_cells mesh
       adjacency */

    cs_real_t  *_xtra_smat = xtra_smat + 2*f_id;
    _xtra_smat[0] = _xtra_smat[1] = contrib;

    /* Diagonal contributions */

    cs_lnum_t cell_i = i_face_cells[f_id][0];
    cs_lnum_t cell_j = i_face_cells[f_id][1];

    diag_smat[cell_i] -= contrib;
    diag_smat[cell_j] -= contrib;

  } /* Loop on interior faces */

  /* Add diagonal contributions from border faces*/

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

    cs_nvec3_t  nvf;
    cs_nvec3(quant->b_face_normal + 3*f_id, &nvf);

    double  contrib = 0;
    for (int k = 0; k < 3; k++)
      contrib += nvf.unitv[k]*nvf.unitv[k];
    contrib *= nvf.meas*nvf.meas;

    /* Diagonal contributions */

    diag_smat[b_face_cells[f_id]] += contrib;

  } /* Loop on border faces */

  /* One assumes a non-symmetric matrix even if in most (all?) cases the matrix
     should be symmetric */

  if (nsp->sles_param->schur_sles_param->solver_class ==
      CS_PARAM_SLES_CLASS_HYPRE)
    sbp->schur_matrix = cs_matrix_external("HYPRE_ParCSR",
                                           false, /* symmetry */
                                           1, 1);
  else
    sbp->schur_matrix = cs_matrix_msr(false, /* symmetry */
                                      1, 1);

  cs_matrix_set_coefficients(sbp->schur_matrix, false, /* symmetry */
                             1, 1,
                             n_i_faces, i_face_cells,
                             diag_smat, xtra_smat);

  /* Return arrays (to be freed when the algorithm is converged) */

  sbp->schur_diag = diag_smat;
  sbp->schur_xtra = xtra_smat;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the matrix for an approximation of the Schur complement based
 *        on the inverse of the sum of the absolute values of the velocity
 *        block
 *
 * \param[in]      nsp     pointer to a cs_navsto_param_t structure
 * \param[in]      ssys    pointer to a saddle-point system structure
 * \param[in, out] sbp     pointer to a cs_saddle_block_precond_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_invlumped_schur_sbp(const cs_navsto_param_t       *nsp,
                     const cs_saddle_system_t      *ssys,
                     cs_saddle_block_precond_t     *sbp)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_lnum_t  n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t  n_i_faces = m->n_i_faces;
  const cs_lnum_t  n_b_faces = m->n_b_faces;
  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  const cs_lnum_t  b11_size = ssys->x1_size;

  /* Compute m11^-1 lumped */

  /* Modify the tolerance in order to be less accurate on this step */

  cs_param_sles_t  *slesp0 = cs_param_sles_create(-1, "schur:inv_lumped");

  cs_param_sles_copy_from(sbp->m11_slesp, slesp0);
  slesp0->eps = 1e-3;
  slesp0->n_max_iter = 50;

  cs_real_t  *rhs = NULL;
  BFT_MALLOC(rhs, b11_size, cs_real_t);
  for (cs_lnum_t i = 0; i < b11_size; i++) rhs[i] = 1;

  cs_real_t  *inv_lumped = NULL;
  BFT_MALLOC(inv_lumped, b11_size, cs_real_t);
  memset(inv_lumped, 0, sizeof(cs_real_t)*b11_size);

  cs_cdo_solve_scalar_system(b11_size,
                             slesp0,
                             ssys->m11_matrices[0],
                             ssys->rset,
                             1,     /* no normalization */
                             false, /* rhs_redux --> already done */
                             sbp->m11_sles,
                             inv_lumped,
                             rhs);

  /* Partial memory free */

  BFT_FREE(rhs);
  cs_param_sles_free(&slesp0);

  /* Native format for the Schur approximation matrix */

  cs_real_t   *diag_smat = NULL;
  cs_real_t   *xtra_smat = NULL;

  BFT_MALLOC(diag_smat, n_cells_ext, cs_real_t);
  BFT_MALLOC(xtra_smat, 2*n_i_faces, cs_real_t);

  memset(diag_smat, 0, n_cells_ext*sizeof(cs_real_t));
  memset(xtra_smat, 0, 2*n_i_faces*sizeof(cs_real_t));

  /* Add diagonal and extra-diagonal contributions from interior faces */

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {

    const cs_real_t  *a_ff = inv_lumped + 3*f_id;
    const cs_nvec3_t  nvf = cs_quant_set_face_nvec(f_id, quant);

    double  contrib = 0;
    for (int k = 0; k < 3; k++)
      contrib += a_ff[k]*nvf.unitv[k]*nvf.unitv[k];
    contrib *= -nvf.meas*nvf.meas;

    /* Extra-diagonal contribution. This is scanned by the i_face_cells mesh
       adjacency */

    cs_real_t  *_xtra_smat = xtra_smat + 2*f_id;
    _xtra_smat[0] = _xtra_smat[1] = contrib;

    /* Diagonal contributions */

    cs_lnum_t cell_i = i_face_cells[f_id][0];
    cs_lnum_t cell_j = i_face_cells[f_id][1];

    diag_smat[cell_i] -= contrib;
    diag_smat[cell_j] -= contrib;

  } /* Loop on interior faces */

  /* Add diagonal contributions from border faces*/

  cs_real_t  *diag_shift = inv_lumped + 3*n_i_faces;
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

    const cs_real_t  *a_ff = diag_shift + 3*f_id;

    cs_nvec3_t  nvf;
    cs_nvec3(quant->b_face_normal + 3*f_id, &nvf);

    double  contrib = 0;
    for (int k = 0; k < 3; k++)
      contrib += a_ff[k]*nvf.unitv[k]*nvf.unitv[k];
    contrib *= nvf.meas*nvf.meas;

    /* Diagonal contributions */

    diag_smat[b_face_cells[f_id]] += contrib;

  } /* Loop on border faces */

  /* Return the associated matrix */

  /* One assumes a non-symmetric matrix even if in most (all?) cases the matrix
     should be symmetric */

  if (nsp->sles_param->schur_sles_param->solver_class ==
      CS_PARAM_SLES_CLASS_HYPRE)
    sbp->schur_matrix = cs_matrix_external("HYPRE_ParCSR",
                                           false, /* symmetry */
                                           1, 1);
  else
    sbp->schur_matrix = cs_matrix_msr(false, /* symmetry */
                                      1, 1);

  cs_matrix_set_coefficients(sbp->schur_matrix, false, /* symmetry */
                             1, 1,
                             n_i_faces, i_face_cells,
                             diag_smat, xtra_smat);

  /* Return arrays (to be freed when the algorithm is converged) */

  sbp->schur_diag = diag_smat;
  sbp->schur_xtra = xtra_smat;

  sbp->m11_inv_diag = inv_lumped;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the matrix for an approximation of the Schur complement based
 *        on the inverse of the sum of the absolute values of the velocity
 *        block
 *
 * \param[in]      nsp         pointer to a cs_navsto_param_t structure
 * \param[in]      ssys        pointer to a saddle-point system structure
 * \param[in]      schur_sles  sles structure dedicated to the Schur complement
 * \param[in, out] sbp         pointer to a cs_saddle_block_precond_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_schur_approximation(const cs_navsto_param_t       *nsp,
                     const cs_saddle_system_t      *ssys,
                     cs_sles_t                     *schur_sles,
                     cs_saddle_block_precond_t     *sbp)
{
  const cs_navsto_param_sles_t  *nslesp = nsp->sles_param;

  cs_param_sles_t  *schur_slesp = nslesp->schur_sles_param;

  sbp->schur_slesp = schur_slesp;
  if (schur_sles == NULL)
    /* This sles structure should have been defined by name */
    sbp->schur_sles = cs_sles_find_or_add(-1, schur_slesp->name);

  /* Compute the schur approximation matrix */

  switch (nslesp->schur_approximation) {

  case CS_PARAM_SCHUR_DIAG_INVERSE:
    _diag_schur_sbp(nsp, ssys, sbp);
    break;
  case CS_PARAM_SCHUR_ELMAN:
    _elman_schur_sbp(nsp, ssys, sbp);
    break;
  case CS_PARAM_SCHUR_IDENTITY:
    break; /* Nothing to do */
  case CS_PARAM_SCHUR_LUMPED_INVERSE:
    _invlumped_schur_sbp(nsp, ssys, sbp);
    break;
  case CS_PARAM_SCHUR_MASS_SCALED:
    _scaled_mass_sbp(nsp, ssys, sbp);
    break; /* Nothing to do */
  case CS_PARAM_SCHUR_MASS_SCALED_DIAG_INVERSE:
    _scaled_mass_sbp(nsp, ssys, sbp);
    _diag_schur_sbp(nsp, ssys, sbp);
    break;
  case CS_PARAM_SCHUR_MASS_SCALED_LUMPED_INVERSE:
    _scaled_mass_sbp(nsp, ssys, sbp);
    _invlumped_schur_sbp(nsp, ssys, sbp);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid Schur approximation.", __func__);
  }
}

#if defined(HAVE_PETSC)
/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup the main iterative solver for the velocity block
 *
 * \param[in]      model    type of model related to the Navsto system
 * \param[in]      nslesp   set of parameter for the monolithic SLES
 * \param[in]      slesp    set of parameters for the velocity SLES
 * \param[in, out] ksp      pointer to PETSc KSP context
 */
/*----------------------------------------------------------------------------*/

static void
_set_petsc_main_solver(const cs_navsto_param_model_t   model,
                       const cs_navsto_param_sles_t   *nslesp,
                       const cs_param_sles_t          *slesp,
                       KSP                             ksp)
{
  if (model == CS_NAVSTO_MODEL_STOKES)
    KSPSetType(ksp, KSPFCG);

  else { /* Advection is present, so one needs a more genric iterative solver */

    /* Flexible GMRES */

    KSPSetType(ksp, KSPFGMRES);
    KSPGMRESSetRestart(ksp, nslesp->il_algo_restart);

  }

  /* Set KSP tolerances */

  PetscReal rtol, abstol, dtol;
  PetscInt  max_it;
  KSPGetTolerances(ksp, &rtol, &abstol, &dtol, &max_it);
  KSPSetTolerances(ksp,
                   nslesp->il_algo_param.rtol,  /* relative convergence tol. */
                   nslesp->il_algo_param.atol,  /* absolute convergence tol. */
                   nslesp->il_algo_param.dtol,  /* divergence tol. */
                   nslesp->il_algo_param.n_max_algo_iter); /* max number iter */

  /* Set the normalization of the residual */

  _set_residual_normalization(slesp->resnorm_type, ksp);
}

#if defined(PETSC_HAVE_HYPRE)
/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup advanced parameters for the AMG related to the velocity field
 *         when BoomerAMG from the HYPRE library is used
 */
/*----------------------------------------------------------------------------*/

static void
_setup_velocity_boomeramg(void)
{
#if PETSC_VERSION_GE(3,7,0)
  PetscOptionsSetValue(NULL,
                       "-pc_velocity_hypre_boomeramg_coarsen_type", "HMIS");
  PetscOptionsSetValue(NULL,
                       "-pc_velocity_hypre_boomeramg_interp_type", "ext+i-cc");
  PetscOptionsSetValue(NULL,
                       "-pc_velocity_hypre_boomeramg_agg_nl", "2");
  PetscOptionsSetValue(NULL,
                       "-pc_velocity_hypre_boomeramg_P_max", "4");
  PetscOptionsSetValue(NULL,
                       "-pc_velocity_hypre_boomeramg_strong_threshold", "0.5");
  PetscOptionsSetValue(NULL,
                       "-pc_velocity_hypre_boomeramg_no_CF", "");
#else
  PetscOptionsSetValue("-pc_velocity_hypre_boomeramg_coarsen_type","HMIS");
  PetscOptionsSetValue("-pc_velocity_hypre_boomeramg_interp_type","ext+i-cc");
  PetscOptionsSetValue("-pc_velocity_hypre_boomeramg_agg_nl","2");
  PetscOptionsSetValue("-pc_velocity_hypre_boomeramg_P_max","4");
  PetscOptionsSetValue("-pc_velocity_hypre_boomeramg_strong_threshold","0.5");
  PetscOptionsSetValue("-pc_velocity_hypre_boomeramg_no_CF","");
#endif
}
#endif  /* PETSC_HAVE_HYPRE */

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup advanced parameters for the AMG related to the velocity field
 *         when GAMG from the PETSc library is used
 */
/*----------------------------------------------------------------------------*/

static void
_setup_velocity_gamg(void)
{
#if PETSC_VERSION_GE(3,7,0)
  PetscOptionsSetValue(NULL, "-mg_velocity_levels_ksp_type", "richardson");
  PetscOptionsSetValue(NULL, "-mg_velocity_levels_pc_type", "sor");
  PetscOptionsSetValue(NULL, "-mg_velocity_levels_ksp_max_it", "1");
  PetscOptionsSetValue(NULL, "-pc_velocity_gamg_threshold", "0.02");
  PetscOptionsSetValue(NULL, "-pc_velocity_gamg_reuse_interpolation", "TRUE");
  PetscOptionsSetValue(NULL, "-pc_velocity_gamg_square_graph", "4");
#else
  PetscOptionsSetValue("-mg_velocity_levels_ksp_type", "richardson");
  PetscOptionsSetValue("-mg_velocity_levels_pc_type", "sor");
  PetscOptionsSetValue("-mg_velocity_levels_ksp_max_it", "1");
  PetscOptionsSetValue("-pc_velocity_gamg_threshold", "0.02");
  PetscOptionsSetValue("-pc_velocity_gamg_reuse_interpolation", "TRUE");
  PetscOptionsSetValue("-pc_velocity_gamg_square_graph", "4");
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Generate IndexSet for the PETSc FieldSplit preconditioner
 *
 * \param[in]       rset    pointer to a range set structure
 * \param[in, out]  isp     IndexSet for the pressure DoFs
 * \param[in, out]  isv     IndexSet for the velocity DoFs
 */
/*----------------------------------------------------------------------------*/

static void
_build_is_for_fieldsplit(const cs_range_set_t   *rset,
                         IS                     *isp,
                         IS                     *isv)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;

  PetscInt  n_faces = quant->n_faces;
  PetscInt  n_cells = quant->n_cells;
  PetscInt  *indices = NULL;

  PetscMalloc1(3*n_faces, &indices);

  /* IndexSet for the velocity DoFs */

  if (rset->n_elts[0] == rset->n_elts[1]) {

    for (PetscInt i = 0; i < 3*n_faces; i++)
      indices[i] = rset->g_id[i];
    ISCreateGeneral(PETSC_COMM_SELF, 3*n_faces, indices, PETSC_COPY_VALUES,
                    isv);

  }
  else {

    PetscInt  n_velocity_elts = 0;
    for (PetscInt i = 0; i < 3*n_faces; i++) {
      cs_gnum_t  g_id = rset->g_id[i];
      if (g_id >= rset->l_range[0] && g_id < rset->l_range[1])
        indices[n_velocity_elts++] = g_id;
    }
    ISCreateGeneral(PETSC_COMM_WORLD, n_velocity_elts, indices,
                    PETSC_COPY_VALUES, isv);

  }

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_SLES_DBG > 1
  /* Print the index set to stdout */
  ISView(*isv, PETSC_VIEWER_STDOUT_SELF);
#endif

  /* Re-used the buffer indices to create the IndexSet for pressure DoFs
   * Pressure unknowns are located at cell centers so the treatment should be
   * the same in sequential and parallel computation */

  for (PetscInt i = 0; i < n_cells; i++)
    indices[i] = rset->g_id[i + 3*n_faces];
  ISCreateGeneral(PETSC_COMM_SELF, n_cells, indices, PETSC_COPY_VALUES,
                  isp);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_SLES_DBG > 1
  /* Print the index set to stdout */
  ISView(*isp, PETSC_VIEWER_STDOUT_SELF);
#endif

  PetscFree(indices);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set command line options for PC according to the kind of
 *        preconditionner
 *
 * \param[in]   slesp      set of parameters for the linear algebra
 */
/*----------------------------------------------------------------------------*/

static PCType
_petsc_get_pc_type(const cs_param_sles_t    *slesp)
{
  PCType  pc_type = PCNONE;

  switch (slesp->precond) {

  case CS_PARAM_PRECOND_NONE:
    return PCNONE;

  case CS_PARAM_PRECOND_DIAG:
    return PCJACOBI;

  case CS_PARAM_PRECOND_BJACOB_ILU0:
  case CS_PARAM_PRECOND_BJACOB_SGS:
    return PCBJACOBI;

  case CS_PARAM_PRECOND_SSOR:
    return PCSOR;

  case CS_PARAM_PRECOND_ICC0:
    return PCICC;

  case CS_PARAM_PRECOND_ILU0:
    return PCILU;

  case CS_PARAM_PRECOND_LU:
    return PCLU;

  case CS_PARAM_PRECOND_AMG:
    {
      switch (slesp->amg_type) {

      case CS_PARAM_AMG_PETSC_GAMG_V:
      case CS_PARAM_AMG_PETSC_GAMG_W:
        return PCGAMG;
        break;

      case CS_PARAM_AMG_PETSC_PCMG:
        return PCMG;
        break;

      case CS_PARAM_AMG_HYPRE_BOOMER_V:
      case CS_PARAM_AMG_HYPRE_BOOMER_W:
#if defined(PETSC_HAVE_HYPRE)
        return PCHYPRE;
#else
        cs_base_warn(__FILE__, __LINE__);
        cs_log_printf(CS_LOG_DEFAULT,
                      "%s: Switch to MG since BoomerAMG is not available.\n",
                      __func__);
        return PCMG;
#endif

      default:
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid AMG type for the PETSc library.", __func__);
        break;

      } /* End of switch on the AMG type */

    } /* AMG as preconditioner */
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Preconditioner not interfaced with PETSc.", __func__);
    break;
  }

  return pc_type;
}

/*----------------------------------------------------------------------------
 * \brief  Function pointer: setup hook for setting PETSc solver and
 *         preconditioner.
 *         Case of additive block preconditioner for a GMRES
 *
 * \param[in]      slesp     pointer to a set of SLES settings
 * \param[in]      rtol      relative tolerance to set
 * \param[in]      max_it    max number of iterations
 * \param[in, out] u_ksp     pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

static void
_set_velocity_ksp(const cs_param_sles_t   *slesp,
                  PetscReal                rtol,
                  PetscInt                 max_it,
                  KSP                      u_ksp)
{
  PC u_pc;
  KSPGetPC(u_ksp, &u_pc);
  PCType  pc_type = _petsc_get_pc_type(slesp);

  _set_residual_normalization(slesp->resnorm_type, u_ksp);

  /* Set the solver */

  switch (slesp->solver) {

  case CS_PARAM_ITSOL_NONE:
    KSPSetType(u_ksp, KSPPREONLY);
    break;
  case CS_PARAM_ITSOL_FCG:
    KSPSetType(u_ksp, KSPFCG);
    break;
  case CS_PARAM_ITSOL_CG:
    KSPSetType(u_ksp, KSPCG);
    break;
  case CS_PARAM_ITSOL_BICG:      /* Improved Bi-CG stab */
    KSPSetType(u_ksp, KSPIBCGS);
    break;
  case CS_PARAM_ITSOL_BICGSTAB2: /* BiCGstab2 */
    KSPSetType(u_ksp, KSPBCGSL);
    break;
  case CS_PARAM_ITSOL_MUMPS:     /* Direct solver (factorization) */
#if defined(PETSC_HAVE_MUMPS)
    KSPSetType(u_ksp, KSPPREONLY);
    PCSetType(u_pc, PCLU);
    PCFactorSetMatSolverType(u_pc, MATSOLVERMUMPS);
#else
    bft_error(__FILE__, __LINE__, 0,
              " %s: MUMPS not interfaced with this installation of PETSc.",
              __func__);
#endif
    break;
  case CS_PARAM_ITSOL_GMRES:
    /* Number of iterations before restarting = 30 (default value)  */
    KSPSetType(u_ksp, KSPGMRES);
    break;
  case CS_PARAM_ITSOL_FGMRES:
    /* Number of iterations before restarting = 30 (default value)  */
    KSPSetType(u_ksp, KSPFGMRES);
    break;

  case CS_PARAM_ITSOL_MUMPS_LDLT:     /* Direct solver (factorization) */
  case CS_PARAM_ITSOL_MUMPS_FLOAT:
  case CS_PARAM_ITSOL_MUMPS_FLOAT_LDLT:
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid solver. Try mumps.",
              __func__);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid solver.", __func__);
    break;

  } /* Switch on solver */

  if (slesp->solver != CS_PARAM_ITSOL_MUMPS)
    PCSetType(u_pc, pc_type);

  /* Additional settings for the preconditioner */

  switch (slesp->precond) {

  case CS_PARAM_PRECOND_AMG:
    switch (slesp->amg_type) {

    case CS_PARAM_AMG_PETSC_GAMG_V:
      PCGAMGSetType(u_pc, PCGAMGAGG);
      PCGAMGSetNSmooths(u_pc, 1);
      PCMGSetCycleType(u_pc, PC_MG_CYCLE_V);
      _setup_velocity_gamg();
      break;

    case CS_PARAM_AMG_PETSC_GAMG_W:
      PCGAMGSetType(u_pc, PCGAMGAGG);
      PCGAMGSetNSmooths(u_pc, 1);
      PCMGSetCycleType(u_pc, PC_MG_CYCLE_W);
      _setup_velocity_gamg();
      break;

    case CS_PARAM_AMG_HYPRE_BOOMER_V:
#if defined(PETSC_HAVE_HYPRE)
      PCHYPRESetType(u_pc, "boomeramg");
      PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_cycle_type","V");
      _setup_velocity_boomeramg();
#else
      PCGAMGSetType(u_pc, PCGAMGAGG);
      PCGAMGSetNSmooths(u_pc, 1);
      PCMGSetCycleType(u_pc, PC_MG_CYCLE_V);
      _setup_velocity_gamg();
#endif
      break;

    case CS_PARAM_AMG_HYPRE_BOOMER_W:
#if defined(PETSC_HAVE_HYPRE)
      PCHYPRESetType(u_pc, "boomeramg");
      PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_cycle_type","W");
      _setup_velocity_boomeramg();
#else
      PCGAMGSetType(u_pc, PCGAMGAGG);
      PCGAMGSetNSmooths(u_pc, 1);
      PCMGSetCycleType(u_pc, PC_MG_CYCLE_W);
      _setup_velocity_gamg();
#endif
      break;

    default:
      bft_error(__FILE__, __LINE__, 0, "%s: Invalid AMG type.", __func__);
      break;

    } /* AMG type */
    break;

  default:
    break; /* Nothing else to do */

  } /* Switch on preconditioner */

  /* Set tolerance and number of iterations */

  PetscReal _rtol, abstol, dtol;
  PetscInt  _max_it;
  KSPGetTolerances(u_ksp, &_rtol, &abstol, &dtol, &_max_it);
  KSPSetTolerances(u_ksp,
                   rtol,        /* relative convergence tolerance */
                   abstol,      /* absolute convergence tolerance */
                   dtol,        /* divergence tolerance */
                   max_it);     /* max number of iterations */

  PCSetFromOptions(u_pc);
  PCSetUp(u_pc);

  KSPSetFromOptions(u_ksp);
  KSPSetUp(u_ksp);
}

/*----------------------------------------------------------------------------
 * \brief  Function pointer: setup hook for setting PETSc solver and
 *         preconditioner.
 *         Case of additive block preconditioner
 *
 * \param[in, out] context     pointer to optional (untyped) value or structure
 * \param[in, out] ksp_struct  pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

static void
_additive_amg_hook(void     *context,
                   void     *ksp_struct)
{
  IS  isv = NULL, isp = NULL;

  KSP  ksp = ksp_struct;
  cs_cdofb_monolithic_petsc_context_t  *phc = context;

  const cs_cdofb_monolithic_t  *sc = phc->sc;
  const cs_navsto_monolithic_t *cc = sc->coupling_context;
  const cs_equation_t  *mom_eq = cc->momentum;

  cs_param_sles_t  *slesp = mom_eq->param->sles_param;

  const cs_navsto_param_t  *nsp = phc->nsp;
  const cs_navsto_param_sles_t  *nslesp = nsp->sles_param;

  cs_fp_exception_disable_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */

  /* Set the main iterative solver for the saddle-point problem */

  _set_petsc_main_solver(nsp->model, nslesp, slesp, ksp);

  /* Apply modifications to the KSP structure */

  PC up_pc, p_pc;

  KSPGetPC(ksp, &up_pc);
  PCSetType(up_pc, PCFIELDSPLIT);
  PCFieldSplitSetType(up_pc, PC_COMPOSITE_ADDITIVE);

  const cs_cdo_system_helper_t  *sh = sc->system_helper;
  assert(sh->n_blocks == 1); /* one block mixing velocity/pressure */

  _build_is_for_fieldsplit(cs_cdo_system_get_range_set(sh, 0),
                           &isp, &isv);

  /* First level Pressure | Velocity (X,Y,Z) */

  PCFieldSplitSetIS(up_pc, "velocity", isv);
  PCFieldSplitSetIS(up_pc, "pressure", isp);

  /* Need to call PCSetUp before configuring the second level (Thanks to
   * Natacha Bereux) */

  PCSetFromOptions(up_pc);
  PCSetUp(up_pc);
  KSPSetUp(ksp);

  PetscInt  n_split;
  KSP  *up_subksp;
  PCFieldSplitGetSubKSP(up_pc, &n_split, &up_subksp);
  assert(n_split == 2);

  KSP  p_ksp = up_subksp[1];
  KSPSetType(p_ksp, KSPPREONLY);
  KSPGetPC(p_ksp, &p_pc);
  PCSetType(p_pc, PCJACOBI);

  PCSetFromOptions(p_pc);
  PCSetUp(p_pc);
  KSPSetUp(p_ksp);

  /* Set the KSP used as preconditioner for the velocity block */

  _set_velocity_ksp(slesp, slesp->eps, slesp->n_max_iter, up_subksp[0]);

  /* User function for additional settings */

  cs_user_sles_petsc_hook(context, ksp);

  /* Apply modifications to the KSP structure */

  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);

  /* Dump the setup related to PETSc in a specific file */

  if (!slesp->setup_done) {
    cs_sles_petsc_log_setup(ksp);
    slesp->setup_done = true;
  }

  PetscFree(up_subksp);
  ISDestroy(&isp);
  ISDestroy(&isv);

  cs_fp_exception_restore_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */
}

/*----------------------------------------------------------------------------
 * \brief  Function pointer: setup hook for setting PETSc solver and
 *         preconditioner.
 *         Case of multiplicative block preconditioner
 *
 * \param[in, out] context    pointer to optional (untyped) value or structure
 * \param[in, out] ksp_struct pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

static void
_multiplicative_hook(void     *context,
                     void     *ksp_struct)
{
  IS  isv = NULL, isp = NULL;

  KSP  ksp = ksp_struct;
  cs_cdofb_monolithic_petsc_context_t  *phc = context;

  const cs_navsto_param_t  *nsp = phc->nsp;
  const cs_navsto_param_sles_t  *nslesp = nsp->sles_param;
  const cs_cdofb_monolithic_t  *sc = phc->sc;
  const cs_navsto_monolithic_t *cc = sc->coupling_context;
  const cs_equation_t  *mom_eq = cc->momentum;

  cs_param_sles_t  *slesp = mom_eq->param->sles_param;

  /* Set the main iterative solver for the saddle-point problem */

  _set_petsc_main_solver(nsp->model, nslesp, slesp, ksp);

  /* Apply modifications to the KSP structure */

  PC up_pc, p_pc;

  KSPGetPC(ksp, &up_pc);
  PCSetType(up_pc, PCFIELDSPLIT);
  PCFieldSplitSetType(up_pc, PC_COMPOSITE_MULTIPLICATIVE);

  const cs_cdo_system_helper_t  *sh = sc->system_helper;
  assert(sh->n_blocks == 1); /* one block mixing velocity/pressure */

  _build_is_for_fieldsplit(cs_cdo_system_get_range_set(sh, 0),
                           &isp, &isv);

  /* First level Pressure | Velocity (X,Y,Z) */

  PCFieldSplitSetIS(up_pc, "velocity", isv);
  PCFieldSplitSetIS(up_pc, "pressure", isp);

  /* Need to call PCSetUp before configuring the second level (Thanks to
     Natacha Bereux) */

  PCSetFromOptions(up_pc);
  PCSetUp(up_pc);
  KSPSetUp(ksp);

  PetscInt  n_split;
  KSP  *up_subksp;
  PCFieldSplitGetSubKSP(up_pc, &n_split, &up_subksp);
  assert(n_split == 2);

  KSP  p_ksp = up_subksp[1];
  KSPSetType(p_ksp, KSPPREONLY);
  KSPGetPC(p_ksp, &p_pc);
  PCSetType(p_pc, PCJACOBI);

  PCSetFromOptions(p_pc);
  PCSetUp(p_pc);
  KSPSetUp(p_ksp);

  /* Set the velocity block */

  _set_velocity_ksp(slesp, slesp->eps, slesp->n_max_iter, up_subksp[0]);

  /* User function for additional settings */

  cs_user_sles_petsc_hook(context, ksp);

  /* Apply modifications to the KSP structure */

  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);

  /* Dump the setup related to PETSc in a specific file */

  if (!slesp->setup_done) {
    cs_sles_petsc_log_setup(ksp);
    slesp->setup_done = true;
  }

  PetscFree(up_subksp);
  ISDestroy(&isp);
  ISDestroy(&isv);

  cs_fp_exception_restore_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */
}

/*----------------------------------------------------------------------------
 * \brief  Function pointer: setup hook for setting PETSc solver and
 *         preconditioner.
 *         Case of diagonal Schur preconditioner by block
 *
 * \param[in, out] context    pointer to optional (untyped) value or structure
 * \param[in, out] ksp_struct pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

static void
_diag_schur_hook(void     *context,
                 void     *ksp_struct)
{
  IS  isv = NULL, isp = NULL;

  KSP  ksp = ksp_struct;
  cs_cdofb_monolithic_petsc_context_t  *phc = context;

  const cs_navsto_param_t  *nsp = phc->nsp;
  const cs_navsto_param_sles_t  *nslesp = nsp->sles_param;
  const cs_cdofb_monolithic_t  *sc = phc->sc;
  const cs_navsto_monolithic_t *cc = sc->coupling_context;
  const cs_equation_t  *mom_eq = cc->momentum;

  cs_param_sles_t  *slesp = mom_eq->param->sles_param;

  cs_fp_exception_disable_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */

  /* Set the main iterative solver for the saddle-point problem */

  _set_petsc_main_solver(nsp->model, nslesp, slesp, ksp);

  /* Apply modifications to the KSP structure */

  PC up_pc, p_pc;

  KSPGetPC(ksp, &up_pc);
  PCSetType(up_pc, PCFIELDSPLIT);
  PCFieldSplitSetType(up_pc, PC_COMPOSITE_SCHUR);
  PCFieldSplitSetSchurFactType(up_pc, PC_FIELDSPLIT_SCHUR_FACT_DIAG);
  PCFieldSplitSetSchurPre(up_pc, PC_FIELDSPLIT_SCHUR_PRE_SELFP, NULL);

  /* Retrieve the range set associated to the block */

  const cs_cdo_system_helper_t  *sh = sc->system_helper;
  assert(sh->n_blocks == 1); /* one block mixing velocity/pressure */

  _build_is_for_fieldsplit(cs_cdo_system_get_range_set(sh, 0),
                           &isp, &isv);

  /* First level Pressure | Velocity (X,Y,Z) */

  PCFieldSplitSetIS(up_pc, "velocity", isv);
  PCFieldSplitSetIS(up_pc, "pressure", isp);

  /* Need to call PCSetUp before configuring the second level (Thanks to
     Natacha Bereux) */

  PCSetFromOptions(up_pc);
  PCSetUp(up_pc);
  KSPSetUp(ksp);

  PetscInt  n_split;
  KSP  *up_subksp;
  PCFieldSplitGetSubKSP(up_pc, &n_split, &up_subksp);
  assert(n_split == 2);

  KSP  p_ksp = up_subksp[1];
  KSPSetType(p_ksp, KSPMINRES);
  KSPGetPC(p_ksp, &p_pc);
  PCSetType(p_pc, PCNONE);

  PCSetFromOptions(p_pc);
  PCSetUp(p_pc);
  KSPSetUp(p_ksp);

  /* Set the velocity block */

  _set_velocity_ksp(slesp, slesp->eps, slesp->n_max_iter, up_subksp[0]);

  /* User function for additional settings */

  cs_user_sles_petsc_hook(context, ksp);

  /* Apply modifications to the KSP structure */

  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);

  /* Dump the setup related to PETSc in a specific file */

  if (!slesp->setup_done) {
    cs_sles_petsc_log_setup(ksp);
    slesp->setup_done = true;
  }

  PetscFree(up_subksp);
  ISDestroy(&isp);
  ISDestroy(&isv);

  cs_fp_exception_restore_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the solver options and its preconditioner in case of Notay's
 *         strategy. Solver is a flexible GMRES by default.
 *
 * \param[in]      slesp   pointer to a cs_sles_param_t structure
 * \param[in, out] ksp     KSP structure to set
 */
/*----------------------------------------------------------------------------*/

static void
_notay_solver(const cs_param_sles_t         *slesp,
              KSP                            ksp)
{
  PC  up_pc;
  KSPGetPC(ksp, &up_pc);

  PCType  pc_type = _petsc_get_pc_type(slesp);
  PCSetType(up_pc, pc_type);

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  cs_gnum_t  system_size = 3*quant->n_g_faces + quant->n_cells;

  /* Additional settings for the preconditioner */
  switch (slesp->precond) {

  case CS_PARAM_PRECOND_AMG:
    switch (slesp->amg_type) {

    case CS_PARAM_AMG_PETSC_GAMG_V:
    case CS_PARAM_AMG_PETSC_GAMG_W:
      _set_gamg_pc("", system_size, slesp->amg_type,
                   false, 1); /* is_sym, smooth_lvl */
      break;

    case CS_PARAM_AMG_HYPRE_BOOMER_V:
#if defined(PETSC_HAVE_HYPRE)
      PCHYPRESetType(up_pc, "boomeramg");
      PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_cycle_type","V");
#else
      _set_gamg_pc("", system_size, slesp->amg_type,
                   false, 1); /* is_sym, smooth_lvl */
#endif
      break;

    case CS_PARAM_AMG_HYPRE_BOOMER_W:
#if defined(PETSC_HAVE_HYPRE)
      PCHYPRESetType(up_pc, "boomeramg");
      PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_cycle_type","W");
#else
      _set_gamg_pc("", system_size, slesp->amg_type,
                   false, 1); /* is_sym, smooth_lvl */
#endif
      break;

    default:
      bft_error(__FILE__, __LINE__, 0, "%s: Invalid AMG type.", __func__);
      break;

    } /* AMG type */
    break;

  default:
    break; /* Nothing else to do */

  } /* Switch on preconditioner */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the solver options with a 2x2 block preconditioner in case of
 *         Notay's strategy
 *
 * \param[in]      rset    pointer to a range set structure
 * \param[in]      slesp   pointer to a cs_sles_param_t structure
 * \param[in, out] ksp     KSP structure to set
 */
/*----------------------------------------------------------------------------*/

static void
_notay_block_precond(const cs_range_set_t       *rset,
                     const cs_param_sles_t      *slesp,
                     KSP                         ksp)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;

  /* Generate an indexed set for each velocity component and the pressure
     block */

  PetscInt  n_faces = quant->n_faces;
  PetscInt  n_cells = quant->n_cells;
  PetscInt  alloc_size = (3*n_faces > n_cells) ? 3*n_faces : n_cells;

  PetscInt  *indices = NULL;
  PetscMalloc1(alloc_size, &indices);

  IS  isv, isp;
  PC  up_pc;
  KSPGetPC(ksp, &up_pc);
  PCSetType(up_pc, PCFIELDSPLIT);

  _build_is_for_fieldsplit(rset, &isp, &isv);

  /* Define the preconditioner */

  PCFieldSplitSetIS(up_pc, "vel", isv);
  PCFieldSplitSetIS(up_pc, "pr", isp);

  switch (slesp->pcd_block_type) {
  case CS_PARAM_PRECOND_BLOCK_UPPER_TRIANGULAR:
  case CS_PARAM_PRECOND_BLOCK_LOWER_TRIANGULAR:
    _petsc_cmd(false, "", "pc_fieldsplit_type", "multiplicative");
    break;

  case CS_PARAM_PRECOND_BLOCK_SYM_GAUSS_SEIDEL:
    _petsc_cmd(false, "", "pc_fieldsplit_type", "symmetric_multiplicative");
    break;

  case CS_PARAM_PRECOND_BLOCK_DIAG:
  default:
    _petsc_cmd(false, "", "pc_fieldsplit_type", "additive");
    break;
  }

  char prefix[2][32] = { "fieldsplit_vel", "fieldsplit_pr" };
  PetscInt  n_split;
  KSP  *up_subksp;
  PCFieldSplitGetSubKSP(up_pc, &n_split, &up_subksp);
  assert(n_split == 2);

  for (int k = 0; k < 2; k++) {

    _petsc_cmd(true, prefix[k], "ksp_type", "preonly");
    _petsc_cmd(true, prefix[k], "ksp_norm_type", "unpreconditioned");

    bool  is_sym = false;
    int  smooth_lvl = 1;
    cs_gnum_t  system_size = 3*quant->n_g_faces;
    if (k == 1) {
      system_size = quant->n_g_cells;
      is_sym = true;
      smooth_lvl = 0;
    }

    /* Additional settings for the preconditioner */

    switch (slesp->precond) {

    case CS_PARAM_PRECOND_AMG:
      switch (slesp->amg_type) {

      case CS_PARAM_AMG_PETSC_GAMG_V:
      case CS_PARAM_AMG_PETSC_GAMG_W:
        _petsc_cmd(true, prefix[k], "pc_type", "gamg");
        _set_gamg_pc(prefix[k], system_size, slesp->amg_type,
                     is_sym, smooth_lvl);
        break;

      case CS_PARAM_AMG_HYPRE_BOOMER_V:
#if defined(PETSC_HAVE_HYPRE)
        _petsc_cmd(true, prefix[k], "pc_type", "hypre");
        _petsc_cmd(true, prefix[k], "pc_hypre_type", "boomeramg");
        _petsc_cmd(true, prefix[k], "pc_hypre_boomeramg_cycle_type","V");
#else
        _petsc_cmd(true, prefix[k], "pc_type", "gamg");
        _set_gamg_pc(prefix[k], system_size, slesp->amg_type,
                     is_sym, smooth_lvl);
#endif
        break;

      case CS_PARAM_AMG_HYPRE_BOOMER_W:
#if defined(PETSC_HAVE_HYPRE)
        _petsc_cmd(true, prefix[k], "pc_type", "hypre");
        _petsc_cmd(true, prefix[k], "pc_hypre_type", "boomeramg");
        _petsc_cmd(true, prefix[k], "pc_hypre_boomeramg_cycle_type","W");
#else
        _petsc_cmd(true, prefix[k], "pc_type", "gamg");
        _set_gamg_pc(prefix[k], system_size, slesp->amg_type,
                     is_sym, smooth_lvl);
#endif
        break;

      default:
        bft_error(__FILE__, __LINE__, 0, "%s: Invalid AMG type.", __func__);
        break;

      } /* AMG type */
      break;

    default:
      break; /* Nothing else to do */

    } /* Switch on preconditioner */

    KSP  _ksp = up_subksp[k];
    PC  _pc;
    KSPGetPC(_ksp, &_pc);
    PCSetFromOptions(_pc);

  } /* Loop on blocks */

  PCSetFromOptions(up_pc);
  PCSetUp(up_pc);

  /* Free temporary memory */

  PetscFree(up_subksp);
  ISDestroy(&isp);
  ISDestroy(&isv);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the solver options with a 4x4 block preconditioner in case of
 *         Notay's strategy
 *
 * \param[in]      rset    pointer to a range set structure
 * \param[in]      slesp   pointer to a cs_sles_param_t structure
 * \param[in, out] ksp     KSP structure to set
 */
/*----------------------------------------------------------------------------*/

static void
_notay_full_block_precond(const cs_range_set_t          *rset,
                          const cs_param_sles_t         *slesp,
                          KSP                            ksp)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;

  /* Generate an indexed set for each velocity component and the pressure
     block */

  PetscInt  n_faces = quant->n_faces;
  PetscInt  n_cells = quant->n_cells;
  PetscInt  alloc_size = (n_faces > n_cells) ? n_faces : n_cells;

  PetscInt  *indices = NULL;
  PetscMalloc1(alloc_size, &indices);

  IS  is[4] = {NULL, NULL, NULL, NULL};
  PC up_pc;
  KSPGetPC(ksp, &up_pc);
  PCSetType(up_pc, PCFIELDSPLIT);

  /* IndexSet for the each component of the velocity DoFs */

  for (int k = 0; k < 3; k++) {

    const cs_gnum_t  *u_g_ids = rset->g_id + k*n_faces;

    if (rset->n_elts[0] == rset->n_elts[1]) { /* Sequential run */

      for (PetscInt i = 0; i < n_faces; i++)
        indices[i] = u_g_ids[i];

      ISCreateGeneral(PETSC_COMM_SELF, n_faces, indices, PETSC_COPY_VALUES,
                      is + k);

    }
    else { /* Parallel run */

      PetscInt  n_elts = 0;
      for (PetscInt i = 0; i < n_faces; i++) {
        if (u_g_ids[i] >= rset->l_range[0] && u_g_ids[i] < rset->l_range[1])
          indices[n_elts++] = u_g_ids[i];
      }

      ISCreateGeneral(PETSC_COMM_WORLD, n_elts, indices, PETSC_COPY_VALUES,
                      is + k);

    }

  } /* Loop on velocity components */

  /* Re-used the buffer indices to create the IndexSet for pressure DoFs
   * Pressure unknowns are located at cell centers so the treatment should be
   * the same in sequential and parallel computation */

  const cs_gnum_t  *pr_g_ids = rset->g_id + 3*n_faces;
  for (PetscInt i = 0; i < n_cells; i++)
    indices[i] = pr_g_ids[i];

  ISCreateGeneral(PETSC_COMM_SELF, n_cells, indices, PETSC_COPY_VALUES, is + 3);

  PetscFree(indices);

  /* Define the preconditioner */

  PCFieldSplitSetIS(up_pc, "velx", is[0]);
  PCFieldSplitSetIS(up_pc, "vely", is[1]);
  PCFieldSplitSetIS(up_pc, "velz", is[2]);
  PCFieldSplitSetIS(up_pc, "pr", is[3]);

  switch (slesp->pcd_block_type) {
  case CS_PARAM_PRECOND_BLOCK_FULL_UPPER_TRIANGULAR:
  case CS_PARAM_PRECOND_BLOCK_FULL_LOWER_TRIANGULAR:
    _petsc_cmd(false, "", "pc_fieldsplit_type", "multiplicative");
    break;

  case CS_PARAM_PRECOND_BLOCK_FULL_SYM_GAUSS_SEIDEL:
    _petsc_cmd(false, "", "pc_fieldsplit_type", "symmetric_multiplicative");
    break;

  case CS_PARAM_PRECOND_BLOCK_FULL_DIAG:
  default:
    _petsc_cmd(false, "", "pc_fieldsplit_type", "additive");
    break;
  }

  char prefix[4][32] = { "fieldsplit_velx",
                         "fieldsplit_vely",
                         "fieldsplit_velz",
                         "fieldsplit_pr" };

  PetscInt  n_split;
  KSP  *up_subksp;
  PCFieldSplitGetSubKSP(up_pc, &n_split, &up_subksp);
  assert(n_split == 4);

  for (int k = 0; k < 4; k++) {

    _petsc_cmd(true, prefix[k], "ksp_type", "preonly");
    _petsc_cmd(true, prefix[k], "ksp_norm_type", "unpreconditioned");

    bool  is_sym = false;
    int  smooth_lvl = 1;
    cs_gnum_t  system_size = quant->n_g_faces;
    if (k == 3) {
      system_size = quant->n_g_cells;
      is_sym = true;
      smooth_lvl = 0;
    }

    /* Additional settings for the preconditioner */

    switch (slesp->precond) {

    case CS_PARAM_PRECOND_AMG:
      switch (slesp->amg_type) {

      case CS_PARAM_AMG_PETSC_GAMG_V:
      case CS_PARAM_AMG_PETSC_GAMG_W:
        _petsc_cmd(true, prefix[k], "pc_type", "gamg");
        _set_gamg_pc(prefix[k], system_size, slesp->amg_type,
                     is_sym, smooth_lvl);
        break;

      case CS_PARAM_AMG_HYPRE_BOOMER_V:
#if defined(PETSC_HAVE_HYPRE)
        _petsc_cmd(true, prefix[k], "pc_type", "hypre");
        _petsc_cmd(true, prefix[k], "pc_hypre_type", "boomeramg");
        _petsc_cmd(true, prefix[k], "pc_hypre_boomeramg_cycle_type","V");
#else
        _petsc_cmd(true, prefix[k], "pc_type", "gamg");
        _set_gamg_pc(prefix[k], system_size, slesp->amg_type,
                     is_sym, smooth_lvl);
#endif
        break;

      case CS_PARAM_AMG_HYPRE_BOOMER_W:
#if defined(PETSC_HAVE_HYPRE)
        _petsc_cmd(true, prefix[k], "pc_type", "hypre");
        _petsc_cmd(true, prefix[k], "pc_hypre_type", "boomeramg");
        _petsc_cmd(true, prefix[k], "pc_hypre_boomeramg_cycle_type","W");
#else
        _petsc_cmd(true, prefix[k], "pc_type", "gamg");
        _set_gamg_pc(prefix[k], system_size, slesp->amg_type,
                     is_sym, smooth_lvl);
#endif
        break;

      default:
        bft_error(__FILE__, __LINE__, 0, "%s: Invalid AMG type.", __func__);
        break;

      } /* AMG type */
      break;

    default:
      break; /* Nothing else to do */

    } /* Switch on preconditioner */

    KSP  _ksp = up_subksp[k];
    PC  _pc;
    KSPGetPC(_ksp, &_pc);
    PCSetFromOptions(_pc);

  } /* Loop on blocks */

  PCSetFromOptions(up_pc);
  PCSetUp(up_pc);

  /* Free temporary memory */

  PetscFree(up_subksp);
  for (int k = 0; k < 4; k++)
    ISDestroy(&(is[k]));
}

/*----------------------------------------------------------------------------
 * \brief  Function pointer: setup hook for setting PETSc solver and
 *         preconditioner.
 *         Case of Notay's transformation.
 *         This relies on the following article.
 *         "Algebraic multigrid for Stokes equations", Y. Notay (2017)
 *         SIAM J. Sci. Comput., Vol. 39 (5), pp 88-111
 *
 * \param[in, out] context    pointer to optional (untyped) value or structure
 * \param[in, out] ksp_struct pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

static void
_notay_hook(void     *context,
            void     *ksp_struct)
{
  IS  isv = NULL, isp = NULL;

  KSP  ksp = ksp_struct;
  cs_cdofb_monolithic_petsc_context_t  *phc = context;

  const cs_navsto_param_t  *nsp = phc->nsp;
  const cs_cdofb_monolithic_t  *sc = phc->sc;
  const cs_navsto_monolithic_t *cc = sc->coupling_context;
  const cs_equation_t  *mom_eq = cc->momentum;

  cs_param_sles_t  *slesp = mom_eq->param->sles_param;
  cs_navsto_param_sles_t  *nslesp = nsp->sles_param;

  cs_fp_exception_disable_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */

  if (cs_glob_n_ranks > 1)
    cs_log_printf(CS_LOG_DEFAULT,
                  " %s (Eq. %s) Warning: Algo. not tested in parallel.\n",
                  __func__, slesp->name);

  /* Build IndexSet structures to extract block matrices */

  const cs_cdo_system_helper_t  *sh = sc->system_helper;
  assert(sh->n_blocks == 1); /* one block mixing velocity/pressure */
  const cs_range_set_t  *rset = cs_cdo_system_get_range_set(sh, 0);

  _build_is_for_fieldsplit(rset, &isp, &isv);

  Mat Amat, Amat_nest;
  KSPGetOperators(ksp, &Amat, NULL);

  /* Retrieve blocks */

  Mat A00, A01, A10;
  MatCreateSubMatrix(Amat, isv, isv, MAT_INITIAL_MATRIX, &A00);
  MatCreateSubMatrix(Amat, isv, isp, MAT_INITIAL_MATRIX, &A01);
  MatCreateSubMatrix(Amat, isp, isv, MAT_INITIAL_MATRIX, &A10);

  PetscInt n_v, n_p;
  MatGetSize(A01, &n_v, &n_p);

  /* Define diag = inv(diag(A00)) */

  Vec diag;
  VecCreate(PETSC_COMM_WORLD, &diag);
  VecSetSizes(diag, PETSC_DECIDE, n_v);
  VecSetType(diag, VECMPI);
  MatGetDiagonal(A00, diag);
  VecReciprocal(diag);

  const PetscReal  alpha = cs_navsto_param_get_notay_scaling();
  if (fabs(alpha - 1.0) > 0)
    VecScale(diag, alpha);

  /* Computing new blocks for the transformed system */

  PetscScalar one = 1.0;

  /* Compute the 01 block
   * First step: temp00 <- Id - A*inv(D_A) */

  Mat temp00;
  MatConvert(A00, MATSAME, MAT_INITIAL_MATRIX, &temp00);

  MatDiagonalScale(temp00, NULL, diag); /* left scaling = NULL;
                                           right scaling = diag */
  MatScale(temp00, -one);

  Vec ones;
  VecCreate(PETSC_COMM_WORLD, &ones);
  VecSetSizes(ones, PETSC_DECIDE, n_v);
  VecSetType(ones, VECMPI);
  VecSet(ones, one);
  MatDiagonalSet(temp00, ones, ADD_VALUES);

  /* temp01 = temp00*A01 = (Id - A*inv(D_A))*Bt */

  Mat temp01;
  MatMatMult(temp00, A01, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp01);

  /* Compute the A11 block
   * A11 = B*inv(D_A)*Bt */

  Mat temp10, A11;;
  MatConvert(A10, MATSAME, MAT_INITIAL_MATRIX, &temp10);
  MatDiagonalScale(temp10, NULL, diag); /* temp10 <- B*inv(D_A) */
  MatMatMult(temp10, A01, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &A11);

  /* Partial free */

  VecDestroy(&diag);
  VecDestroy(&ones);
  MatDestroy(&A01);
  MatDestroy(&temp00);
  MatDestroy(&temp10);

  /* Compute A10 <- -1.0*B */

  MatScale(A10, -1.0);

  /* Update blocks and assemble Amat */

  Mat subA[4] = {A00, temp01, A10, A11};

  MatCreateNest(PETSC_COMM_WORLD, 2, NULL, 2, NULL, subA, &Amat_nest);
  MatConvert(Amat_nest, MATMPIAIJ, MAT_INITIAL_MATRIX, &Amat);

  KSPSetOperators(ksp, Amat, Amat);

  /* Partial free */

  MatDestroy(&Amat_nest);
  MatDestroy(&A00);
  MatDestroy(&A10);
  MatDestroy(&A11);
  MatDestroy(&temp01);
  ISDestroy(&isp);
  ISDestroy(&isv);
  MatDestroy(&Amat);

  PC  up_pc;
  KSPGetPC(ksp, &up_pc);

  /* Set the main solver and main options for the preconditioner */

  switch (slesp->solver) {

  case CS_PARAM_ITSOL_MUMPS:
  case CS_PARAM_ITSOL_MUMPS_FLOAT:
#if defined(PETSC_HAVE_MUMPS)
    {
      KSPSetType(ksp, KSPPREONLY);
      PCSetType(up_pc, PCLU);
      PCFactorSetMatSolverType(up_pc, MATSOLVERMUMPS);
    }
#else
    bft_error(__FILE__, __LINE__, 0,
              " %s: MUMPS not interfaced with this installation of PETSc.",
              __func__);
#endif
    break;

  default:
    {
      if (slesp->solver != CS_PARAM_ITSOL_FGMRES) {
        cs_base_warn(__FILE__, __LINE__);
        cs_log_printf(CS_LOG_DEFAULT,
                      "%s: Switch to FGMRES solver.\n", __func__);
      }

      /* Set the solver parameters */

      KSPSetType(ksp, KSPFGMRES);
      KSPGMRESSetRestart(ksp, nslesp->il_algo_restart);

      /* Set KSP tolerances */

      PetscReal rtol, abstol, dtol;
      PetscInt  max_it;
      KSPGetTolerances(ksp, &rtol, &abstol, &dtol, &max_it);
      KSPSetTolerances(ksp,
                       nslesp->il_algo_param.rtol,      /* relative tol. */
                       nslesp->il_algo_param.atol,      /* absolute tol. */
                       nslesp->il_algo_param.dtol,      /* divergence tol. */
                       nslesp->il_algo_param.n_max_algo_iter);

      KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED);

      switch (slesp->pcd_block_type) {

      case CS_PARAM_PRECOND_BLOCK_DIAG:
      case CS_PARAM_PRECOND_BLOCK_LOWER_TRIANGULAR:
      case CS_PARAM_PRECOND_BLOCK_SYM_GAUSS_SEIDEL:
      case CS_PARAM_PRECOND_BLOCK_UPPER_TRIANGULAR:
        _notay_block_precond(rset, slesp, ksp);
        break;

      case CS_PARAM_PRECOND_BLOCK_FULL_DIAG:
      case CS_PARAM_PRECOND_BLOCK_FULL_LOWER_TRIANGULAR:
      case CS_PARAM_PRECOND_BLOCK_FULL_SYM_GAUSS_SEIDEL:
      case CS_PARAM_PRECOND_BLOCK_FULL_UPPER_TRIANGULAR:
        _notay_full_block_precond(rset, slesp, ksp);
        break;

      default:
        _notay_solver(slesp, ksp);
        break;

      } /* Switch on type of block preconditionner */

    } /* Default case */
    break;

  } /* End of switch on the main solver */

  /* User function for additional settings */

  cs_user_sles_petsc_hook(context, ksp);

  /* Apply modifications to the KSP structure */

  PCSetFromOptions(up_pc);
  PCSetUp(up_pc);
  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);

  /* Dump the setup related to PETSc in a specific file */

  if (!slesp->setup_done) {
    cs_sles_petsc_log_setup(ksp);
    slesp->setup_done = true;
  }

  cs_fp_exception_restore_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */
}

/*----------------------------------------------------------------------------
 * \brief  Function pointer: setup hook for setting PETSc solver and
 *         preconditioner.
 *         Case of upper Schur preconditioner by block
 *
 * \param[in, out] context    pointer to optional (untyped) value or structure
 * \param[in, out] ksp_struct pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

static void
_upper_schur_hook(void     *context,
                  void     *ksp_struct)
{
  IS  isv = NULL, isp = NULL;

  KSP  ksp = ksp_struct;
  cs_cdofb_monolithic_petsc_context_t  *phc = context;

  const cs_navsto_param_t  *nsp = phc->nsp;
  const cs_cdofb_monolithic_t  *sc = phc->sc;
  const cs_navsto_monolithic_t *cc = sc->coupling_context;
  const cs_equation_t  *mom_eq = cc->momentum;

  cs_param_sles_t  *slesp = mom_eq->param->sles_param;
  cs_navsto_param_sles_t  *nslesp = nsp->sles_param;

  cs_fp_exception_disable_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */

  /* Set the main iterative solver for the saddle-point problem */

  _set_petsc_main_solver(nsp->model, nslesp, slesp, ksp);

  /* Apply modifications to the KSP structure */

  PC up_pc, p_pc;

  KSPGetPC(ksp, &up_pc);
  PCSetType(up_pc, PCFIELDSPLIT);
  PCFieldSplitSetType(up_pc, PC_COMPOSITE_SCHUR);
  PCFieldSplitSetSchurFactType(up_pc, PC_FIELDSPLIT_SCHUR_FACT_UPPER);
  PCFieldSplitSetSchurPre(up_pc, PC_FIELDSPLIT_SCHUR_PRE_SELFP, NULL);

  const cs_cdo_system_helper_t  *sh = sc->system_helper;
  assert(sh->n_blocks == 1); /* one block mixing velocity/pressure */

  _build_is_for_fieldsplit(cs_cdo_system_get_range_set(sh, 0),
                           &isp, &isv);

  /* First level Pressure | Velocity (X,Y,Z) */

  PCFieldSplitSetIS(up_pc, "velocity", isv);
  PCFieldSplitSetIS(up_pc, "pressure", isp);

  /* Need to call PCSetUp before configuring the second level (Thanks to
     Natacha Bereux) */

  PCSetFromOptions(up_pc);
  PCSetUp(up_pc);
  KSPSetUp(ksp);

  PetscInt  n_split;
  KSP  *up_subksp;
  PCFieldSplitGetSubKSP(up_pc, &n_split, &up_subksp);
  assert(n_split == 2);

  KSP  p_ksp = up_subksp[1];
  KSPSetType(p_ksp, KSPMINRES);
  KSPGetPC(p_ksp, &p_pc);
  PCSetType(p_pc, PCNONE);

  PCSetFromOptions(p_pc);
  PCSetUp(p_pc);
  KSPSetUp(p_ksp);

  /* Set the velocity block */

  _set_velocity_ksp(slesp, slesp->eps, slesp->n_max_iter, up_subksp[0]);

  /* User function for additional settings */

  cs_user_sles_petsc_hook(context, ksp);

  /* Apply modifications to the KSP structure */

  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);

  /* Dump the setup related to PETSc in a specific file */

  if (!slesp->setup_done) {
    cs_sles_petsc_log_setup(ksp);
    slesp->setup_done = true;
  }

  PetscFree(up_subksp);
  ISDestroy(&isp);
  ISDestroy(&isv);

  cs_fp_exception_restore_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */
}

#if PETSC_VERSION_GE(3,11,0)
/*----------------------------------------------------------------------------
 * \brief  Function pointer: setup hook for setting PETSc solver and
 *         preconditioner.
 *         Case of GKB as a solver.
 *
 * \param[in, out] context    pointer to optional (untyped) value or structure
 * \param[in, out] ksp_struct pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

static void
_gkb_hook(void     *context,
          void     *ksp_struct)
{
  IS  isv = NULL, isp = NULL;

  KSP  ksp = ksp_struct;
  cs_cdofb_monolithic_petsc_context_t  *phc = context;

  const cs_navsto_param_t  *nsp = phc->nsp;
  const cs_cdofb_monolithic_t  *sc = phc->sc;
  const cs_navsto_monolithic_t *cc = sc->coupling_context;
  const cs_equation_t  *mom_eq = cc->momentum;

  cs_param_sles_t  *slesp = mom_eq->param->sles_param;
  cs_navsto_param_sles_t  *nslesp = nsp->sles_param;

  cs_fp_exception_disable_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */

  KSPSetType(ksp, KSPPREONLY);

  /* Apply modifications to the KSP structure */

  PC up_pc;

  KSPGetPC(ksp, &up_pc);
  PCSetType(up_pc, PCFIELDSPLIT);
  PCFieldSplitSetType(up_pc, PC_COMPOSITE_GKB);

  PCFieldSplitSetGKBTol(up_pc, 10*nslesp->il_algo_param.rtol);
  PCFieldSplitSetGKBMaxit(up_pc, nslesp->il_algo_param.n_max_algo_iter);
  PCFieldSplitSetGKBNu(up_pc, 0);
  PCFieldSplitSetGKBDelay(up_pc, CS_GKB_TRUNCATION_THRESHOLD);

  const cs_cdo_system_helper_t  *sh = sc->system_helper;
  assert(sh->n_blocks == 1); /* one block mixing velocity/pressure */

  _build_is_for_fieldsplit(cs_cdo_system_get_range_set(sh, 0),
                           &isp, &isv);

  /* First level Pressure | Velocity (X,Y,Z) */

  PCFieldSplitSetIS(up_pc, "velocity", isv);
  PCFieldSplitSetIS(up_pc, "pressure", isp);

  /* Need to call PCSetUp before configuring the second level (Thanks to
     Natacha Bereux) */

  PCSetFromOptions(up_pc);
  PCSetUp(up_pc);
  KSPSetUp(ksp);

  PetscInt  n_split;
  KSP  *up_subksp;
  PCFieldSplitGetSubKSP(up_pc, &n_split, &up_subksp);
  assert(n_split == 2);

  _set_velocity_ksp(slesp, slesp->eps, slesp->n_max_iter, up_subksp[0]);

  /* User function for additional settings */

  cs_user_sles_petsc_hook(context, ksp);

  /* Apply modifications to the KSP structure */

  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);

  /* Dump the setup related to PETSc in a specific file */

  if (!slesp->setup_done) {
    cs_sles_petsc_log_setup(ksp);
    slesp->setup_done = true;
  }

  PetscFree(up_subksp);
  ISDestroy(&isp);
  ISDestroy(&isv);

  cs_fp_exception_restore_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */
}

/*----------------------------------------------------------------------------
 * \brief  Function pointer: setup hook for setting PETSc solver and
 *         preconditioner.
 *         Case of GKB as a preconditioner.
 *
 * \param[in, out] context    pointer to optional (untyped) value or structure
 * \param[in, out] ksp_struct pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

static void
_gkb_precond_hook(void     *context,
                  void     *ksp_struct)
{
  IS  isv = NULL, isp = NULL;

  KSP  ksp = ksp_struct;
  cs_cdofb_monolithic_petsc_context_t  *phc = context;

  const cs_navsto_param_t  *nsp = phc->nsp;
  const cs_cdofb_monolithic_t  *sc = phc->sc;
  const cs_navsto_monolithic_t *cc = sc->coupling_context;
  const cs_equation_t  *mom_eq = cc->momentum;

  cs_param_sles_t  *slesp = mom_eq->param->sles_param;
  cs_navsto_param_sles_t  *nslesp = nsp->sles_param;

  cs_fp_exception_disable_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */

  /* Set the main iterative solver for the saddle-point */

  _set_petsc_main_solver(nsp->model, nslesp, slesp, ksp);

  /* Apply modifications to the KSP structure */

  PC up_pc;

  KSPGetPC(ksp, &up_pc);
  PCSetType(up_pc, PCFIELDSPLIT);
  PCFieldSplitSetType(up_pc, PC_COMPOSITE_GKB);

  /* Default settings for the GKB as preconditioner */

  PCFieldSplitSetGKBTol(up_pc,  1e-2);
  PCFieldSplitSetGKBMaxit(up_pc, 10);
  PCFieldSplitSetGKBNu(up_pc, 0); /* No augmentation */
  PCFieldSplitSetGKBDelay(up_pc, CS_GKB_TRUNCATION_THRESHOLD);

  const cs_cdo_system_helper_t  *sh = sc->system_helper;
  assert(sh->n_blocks == 1); /* one block mixing velocity/pressure */

  _build_is_for_fieldsplit(cs_cdo_system_get_range_set(sh, 0),
                           &isp, &isv);

  /* First level Pressure | Velocity (X,Y,Z) */

  PCFieldSplitSetIS(up_pc, "velocity", isv);
  PCFieldSplitSetIS(up_pc, "pressure", isp);

  /* Need to call PCSetUp before configuring the second level (Thanks to
     Natacha Bereux) */

  PCSetFromOptions(up_pc);
  PCSetUp(up_pc);
  KSPSetUp(ksp);

  PetscInt  n_split;
  KSP  *up_subksp;
  PCFieldSplitGetSubKSP(up_pc, &n_split, &up_subksp);
  assert(n_split == 2);

  /* Set KSP options for the velocity block */

  _set_velocity_ksp(slesp, slesp->eps, slesp->n_max_iter, up_subksp[0]);

  /* User function for additional settings */

  cs_user_sles_petsc_hook(context, ksp);

  /* Apply modifications to the KSP structure */

  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);

  /* Dump the setup related to PETSc in a specific file */

  if (!slesp->setup_done) {
    cs_sles_petsc_log_setup(ksp);
    slesp->setup_done = true;
  }

  PetscFree(up_subksp);
  ISDestroy(&isp);
  ISDestroy(&isv);

  cs_fp_exception_restore_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */
}
#endif  /* GKB available only if version >= 3.11 */

#if defined(PETSC_HAVE_MUMPS)
/*----------------------------------------------------------------------------
 * \brief  Function pointer: setup hook for setting PETSc solver and
 *         preconditioner.
 *         Case of MUMPS via PETSc
 *
 * \param[in, out] context  pointer to optional (untyped) value or structure
 * \param[in, out] ksp      pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

static void
_mumps_hook(void     *context,
            KSP       ksp)
{
  cs_equation_param_t  *eqp = (cs_equation_param_t *)context;
  cs_param_sles_t  *slesp = eqp->sles_param;

  cs_fp_exception_disable_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */

  PC  pc;
  KSPSetType(ksp, KSPPREONLY);
  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCLU);
  PCFactorSetMatSolverType(pc, MATSOLVERMUMPS);

  PetscReal rtol, abstol, dtol;
  PetscInt  max_it;
  KSPGetTolerances(ksp, &rtol, &abstol, &dtol, &max_it);
  KSPSetTolerances(ksp,
                   slesp->eps,  /* relative convergence tolerance */
                   abstol,      /* absolute convergence tolerance */
                   dtol,        /* divergence tolerance */
                   slesp->n_max_iter); /* max number of iterations */

  /* User function for additional settings */

  cs_user_sles_petsc_hook(context, ksp);

  /* Dump the setup related to PETSc in a specific file */

  if (!slesp->setup_done) {
    cs_sles_petsc_log_setup(ksp);
    slesp->setup_done = true;
  }

  cs_fp_exception_restore_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */
}
#endif  /* PETSC_HAVE_MUMPS */
#endif  /* HAVE_PETSC */

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and initialize a GKB builder structure
 *
 * \param[in]  nsp        pointer to a cs_navsto_param_t structure
 * \param[in]  gamma      value of the grad-div coefficient
 * \param[in]  n_u_dofs   number of velocity DoFs (degrees of freedom)
 * \param[in]  n_p_dofs   number of pressure DoFs
 *
 * \return a pointer to a new allocated GKB builder
 */
/*----------------------------------------------------------------------------*/

static cs_gkb_builder_t *
_init_gkb_builder(const cs_navsto_param_t    *nsp,
                  cs_real_t                   gamma,
                  cs_lnum_t                   n_u_dofs,
                  cs_lnum_t                   n_p_dofs)
{
  cs_gkb_builder_t  *gkb = NULL;

  BFT_MALLOC(gkb, 1, cs_gkb_builder_t);

  gkb->gamma = gamma;
  gkb->n_u_dofs = n_u_dofs;
  gkb->n_p_dofs = n_p_dofs;

  /* Vector transformation */

  BFT_MALLOC(gkb->u_tilda, n_u_dofs, cs_real_t);

  /* Rk: b_tilda stores quantities in space M and N alternatively */

  assert(n_u_dofs >= n_p_dofs);
  BFT_MALLOC(gkb->b_tilda, n_u_dofs, cs_real_t);

  /* Auxiliary vectors */

  BFT_MALLOC(gkb->v, n_u_dofs, cs_real_t);
  memset(gkb->v, 0, n_u_dofs*sizeof(cs_real_t));

  BFT_MALLOC(gkb->q, n_p_dofs, cs_real_t);
  BFT_MALLOC(gkb->d, n_p_dofs, cs_real_t);
  BFT_MALLOC(gkb->d__v, n_p_dofs, cs_real_t);
  BFT_MALLOC(gkb->dt_q, n_u_dofs, cs_real_t);
  BFT_MALLOC(gkb->m__v, n_u_dofs, cs_real_t);

  /* Orthogonalization coefficients */

  gkb->alpha = gkb->beta = gkb->zeta = 0.;

  /* Convergence members */

  if (gamma < 1)
    gkb->z_size = CS_GKB_TRUNCATION_THRESHOLD + 1;
  else if (gamma < 10)
    gkb->z_size = CS_GKB_TRUNCATION_THRESHOLD;
  else if (gamma < 100)
    gkb->z_size = CS_MAX(1, CS_GKB_TRUNCATION_THRESHOLD - 1);
  else if (gamma < 1e3)
    gkb->z_size = CS_MAX(1, CS_GKB_TRUNCATION_THRESHOLD - 2);
  else if (gamma < 1e4)
    gkb->z_size = CS_MAX(1, CS_GKB_TRUNCATION_THRESHOLD - 3);
  else
    gkb->z_size = CS_MAX(1, CS_GKB_TRUNCATION_THRESHOLD - 4);

  BFT_MALLOC(gkb->zeta_array, gkb->z_size, cs_real_t);
  memset(gkb->zeta_array, 0, gkb->z_size*sizeof(cs_real_t));

  gkb->zeta_square_sum = 0.;

  const cs_navsto_param_sles_t  *nslesp = nsp->sles_param;

  gkb->algo = cs_iter_algo_create(nslesp->il_algo_param);

  return gkb;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a GKB builder structure
 *
 * \param[in, out]  p_gkb   double pointer to a GKB builder structure
 */
/*----------------------------------------------------------------------------*/

static void
_free_gkb_builder(cs_gkb_builder_t   **p_gkb)
{
  cs_gkb_builder_t  *gkb = *p_gkb;

  if (gkb == NULL)
    return;

  BFT_FREE(gkb->b_tilda);
  BFT_FREE(gkb->u_tilda);

  BFT_FREE(gkb->q);
  BFT_FREE(gkb->d);
  BFT_FREE(gkb->d__v);
  BFT_FREE(gkb->dt_q);
  BFT_FREE(gkb->m__v);
  BFT_FREE(gkb->v);

  BFT_FREE(gkb->zeta_array);

  BFT_FREE(gkb->algo);

  BFT_FREE(gkb);
  *p_gkb = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and initialize a Uzawa builder structure
 *
 * \param[in]  nsp        pointer to a cs_navsto_param_t structure
 * \param[in]  gamma      value of the grad-div coefficient
 * \param[in]  n_u_dofs   number of velocity DoFs (degrees of freedom)
 * \param[in]  n_p_dofs   number of pressure DoFs
 * \param[in]  quant      pointer to additional mesh quantities
 *
 * \return a pointer to a new allocated Uzawa builder
 */
/*----------------------------------------------------------------------------*/

static cs_uza_builder_t *
_init_uzawa_builder(const cs_navsto_param_t      *nsp,
                    cs_real_t                     gamma,
                    cs_lnum_t                     n_u_dofs,
                    cs_lnum_t                     n_p_dofs,
                    const cs_cdo_quantities_t    *quant)
{
  cs_uza_builder_t  *uza = NULL;

  BFT_MALLOC(uza, 1, cs_uza_builder_t);

  uza->alpha = 0;
  uza->gamma = gamma;
  uza->n_u_dofs = n_u_dofs;
  uza->n_p_dofs = n_p_dofs;

  BFT_MALLOC(uza->b_tilda, n_u_dofs, cs_real_t);

  /* Auxiliary vectors */

  BFT_MALLOC(uza->inv_mp, n_p_dofs, cs_real_t);
  BFT_MALLOC(uza->res_p, n_p_dofs, cs_real_t);
  BFT_MALLOC(uza->d__v, n_p_dofs, cs_real_t);
  BFT_MALLOC(uza->rhs, n_u_dofs, cs_real_t);

  uza->gk = NULL;
  uza->dzk = NULL;

  if (nsp->sles_param->strategy == CS_NAVSTO_SLES_UZAWA_CG) {

    const cs_time_step_t  *ts = cs_glob_time_step;

    /* Since gk is used as a variable in a cell system, one has to take into
       account extra-space for synchronization */

    cs_lnum_t  size = n_p_dofs;
    if (cs_glob_n_ranks > 1)
      size = CS_MAX(n_p_dofs, cs_glob_mesh->n_cells_with_ghosts);
    BFT_MALLOC(uza->gk, size, cs_real_t);

    BFT_MALLOC(uza->dzk, n_u_dofs, cs_real_t);

    /* Define alpha weighting */

    cs_real_t  alpha = 1.;
    if (nsp->model_flag & CS_NAVSTO_MODEL_STEADY)
      alpha = 0.01*nsp->lam_viscosity->ref_value;
    else
      alpha /= ts->dt[0];

    const cs_real_t  rho0 = nsp->mass_density->ref_value;
    uza->alpha = rho0*alpha;

    /* Define the inverse of the pressure mass matrix scaled by the viscosity */

    cs_real_t  *visc_val = NULL;
    int  visc_stride = 0;

    if (nsp->turbulence->model->iturb == CS_TURB_NONE) {
      BFT_MALLOC(visc_val, 1, cs_real_t);
      visc_val[0] = nsp->lam_viscosity->ref_value;
    }
    else {
      visc_stride = 1;
      BFT_MALLOC(visc_val, n_p_dofs, cs_real_t);
      cs_property_eval_at_cells(ts->t_cur, nsp->tot_viscosity, visc_val);
    }

    for (cs_lnum_t i = 0; i < n_p_dofs; i++)
      uza->inv_mp[i] = visc_val[visc_stride*i]/quant->cell_vol[i];

    BFT_FREE(visc_val);

  }
  else  {

#   pragma omp parallel for if (uza->n_p_dofs > CS_THR_MIN)
    for (cs_lnum_t ip = 0; ip < uza->n_p_dofs; ip++)
      uza->inv_mp[ip] = 1./quant->cell_vol[ip];

  }

  const cs_navsto_param_sles_t  *nslesp = nsp->sles_param;

  uza->algo = cs_iter_algo_create(nslesp->il_algo_param);

  return uza;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a Uzawa builder structure
 *
 * \param[in, out]  p_uza   double pointer to a Uzawa builder structure
 */
/*----------------------------------------------------------------------------*/

static void
_free_uza_builder(cs_uza_builder_t   **p_uza)
{
  cs_uza_builder_t  *uza = *p_uza;

  if (uza == NULL)
    return;

  BFT_FREE(uza->b_tilda);

  BFT_FREE(uza->inv_mp);
  BFT_FREE(uza->res_p);
  BFT_FREE(uza->d__v);
  BFT_FREE(uza->rhs);
  BFT_FREE(uza->gk);
  BFT_FREE(uza->dzk);

  BFT_FREE(uza->algo);

  BFT_FREE(uza);
  *p_uza = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Apply the divergence operator and store the result in div_v
 *
 * \param[in]      div_op  pointer to the values of divergence operator
 * \param[in]      v       vector to apply in velocity space
 * \param[in, out] div_v   resulting vector in pressure space
 */
/*----------------------------------------------------------------------------*/

static void
_apply_div_op(const cs_real_t   *div_op,
              const cs_real_t   *v,
              cs_real_t         *div_v)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_adjacency_t  *c2f = cs_shared_connect->c2f;

# pragma omp parallel for if (quant->n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

    cs_real_t _div_v = 0;
    for (cs_lnum_t j = c2f->idx[c_id]; j < c2f->idx[c_id+1]; j++)
      _div_v += cs_math_3_dot_product(div_op + 3*j, v + 3*c2f->ids[j]);
    div_v[c_id] = _div_v;

  } /* Loop on cells */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Apply the gradient operator (which is the transpose of the
 *         divergence operator) and store the result in dt_q
 *
 * \param[in]      div_op   pointer to the values of divergence operator
 * \param[in]      q        vector to apply in pressure space
 * \param[in, out] dt_q     resulting vector in velocity space
 */
/*----------------------------------------------------------------------------*/

static void
_apply_div_op_transpose(const cs_real_t   *div_op,
                        const cs_real_t   *q,
                        cs_real_t         *dt_q)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_adjacency_t  *c2f = cs_shared_connect->c2f;

  memset(dt_q, 0, 3*quant->n_faces*sizeof(cs_real_t));

# pragma omp parallel for if (quant->n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

    const cs_real_t  qc = q[c_id];
    for (cs_lnum_t j = c2f->idx[c_id]; j < c2f->idx[c_id+1]; j++) {

      const cs_real_t  *_div_f = div_op + 3*j;

      cs_real_t  *_dt_q = dt_q + 3*c2f->ids[j];
#     pragma omp critical
      {
        _dt_q[0] += qc * _div_f[0];
        _dt_q[1] += qc * _div_f[1];
        _dt_q[2] += qc * _div_f[2];
      }

    } /* Loop on cell faces */

  } /* Loop on cells */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Transform the initial saddle-point problem. The velocity unknown
 *         is modified and is stored in u_tilda as well as the RHS related to
 *         the mass equation and stored in b_tilda
 *
 * \param[in]      matrix   pointer to a cs_matrix_t structure
 * \param[in]      rset     pointer to a range set structure
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in]      nslesp   pointer to SLES settings for NavSto
 * \param[in]      div_op   pointer to the values of divergence operator
 * \param[in, out] slesp    pointer to a set of parameters to drive the SLES
 * \param[in, out] gkb      pointer to a GKB builder structure
 * \param[in, out] sles     pointer to a cs_sles_t structure
 * \param[in]      u_f      initial velocity on faces
 * \param[in]      b_f      right-hand side (scatter/gather if needed) on faces
 * \param[in]      b_c      right_hand side on cells (mass equation)
 */
/*----------------------------------------------------------------------------*/

static void
_transform_gkb_system(const cs_matrix_t              *matrix,
                      const cs_range_set_t           *rset,
                      const cs_equation_param_t      *eqp,
                      const cs_navsto_param_sles_t   *nslesp,
                      const cs_real_t                *div_op,
                      cs_param_sles_t                *slesp,
                      cs_gkb_builder_t               *gkb,
                      cs_sles_t                      *sles,
                      const cs_real_t                *u_f,
                      const cs_real_t                *b_f,
                      const cs_real_t                *b_c)
{
  assert(gkb != NULL);

  cs_real_t  normalization = 1.0; /* TODO */

  if (gkb->gamma > 0) {

#   pragma omp parallel for if (gkb->n_p_dofs > CS_THR_MIN)
    for (cs_lnum_t ip = 0; ip < gkb->n_p_dofs; ip++)
      gkb->b_tilda[ip] = gkb->gamma*b_c[ip]/cs_shared_quant->cell_vol[ip];

    /* Build Dt.b_tilda */

    _apply_div_op_transpose(div_op, gkb->b_tilda, gkb->dt_q);

#   pragma omp parallel for if (gkb->n_u_dofs > CS_THR_MIN)
    for (cs_lnum_t iu = 0; iu < gkb->n_u_dofs; iu++)
      gkb->b_tilda[iu] = b_f[iu] + gkb->dt_q[iu];

  }
  else
    memcpy(gkb->b_tilda, b_f, gkb->n_u_dofs*sizeof(cs_real_t));

  /* Modifiy the tolerance in order to be more accurate on the next solve
     step (the final accuracy relies on this step) */

  char  *init_system_name = slesp->name;
  double  init_eps = slesp->eps;

  char  *system_name = NULL;
  BFT_MALLOC(system_name, strlen(eqp->name) + strlen(":gkb_transfo") + 1, char);
  sprintf(system_name, "%s:gkb_transfo", eqp->name);

  slesp->name = system_name;
  slesp->eps = nslesp->il_algo_param.rtol;

  /* Compute M^-1.(b_f + gamma. Bt.N^-1.b_c) */

  gkb->algo->n_inner_iter
    += (gkb->algo->last_inner_iter
        = cs_cdo_solve_scalar_system(gkb->n_u_dofs,
                                     slesp,
                                     matrix,
                                     rset,
                                     normalization,
                                     true, /* rhs_redux, */
                                     sles,
                                     gkb->v,
                                     gkb->b_tilda));

  /* Set back the initial parameters */

  slesp->name = init_system_name;
  slesp->eps = init_eps;

  /* Compute the initial u_tilda := u_f - M^-1.(b_f + gamma. Bt.N^-1.b_c) */

# pragma omp parallel for if (gkb->n_u_dofs > CS_THR_MIN)
  for (cs_lnum_t iu = 0; iu < gkb->n_u_dofs; iu++)
    gkb->u_tilda[iu] = u_f[iu] - gkb->v[iu];

  /* Compute b_tilda := b_c - div(M^-1.b_f) */

  _apply_div_op(div_op, gkb->v, gkb->d__v);

# pragma omp parallel for if (gkb->n_p_dofs > CS_THR_MIN)
  for (cs_lnum_t ip = 0; ip < gkb->n_p_dofs; ip++)
    gkb->b_tilda[ip] = b_c[ip] - gkb->d__v[ip];

  /* Free memory */

  BFT_FREE(system_name);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the GKB algorithm
 *
 * \param[in]      matrix   pointer to a cs_matrix_t structure
 * \param[in]      rset     pointer to a range set structure
 * \param[in]      div_op   pointer to the values of divergence operator
 * \param[in, out] slesp    pointer to a set of parameters to drive the SLES
 * \param[in, out] gkb      pointer to a GKB builder structure
 * \param[in, out] sles     pointer to a cs_sles_t structure
 * \param[in, out] p_c      right_hand side on cells (mass equation)
 */
/*----------------------------------------------------------------------------*/

static void
_init_gkb_algo(const cs_matrix_t             *matrix,
               const cs_range_set_t          *rset,
               const cs_real_t               *div_op,
               cs_param_sles_t               *slesp,
               cs_gkb_builder_t              *gkb,
               cs_sles_t                     *sles,
               cs_real_t                     *p_c)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  size = quant->n_cells;

  double beta2 = 0.0;

  /* Compute beta := ||b_tilta||_N^-1 and q := N^-1(b_tilda)/beta */

# pragma omp parallel reduction(+:beta2) if (size > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    _thread_range(size, &s_id, &e_id);

    const cs_lnum_t  n = e_id - s_id;
    const cs_real_t  *_w = quant->cell_vol + s_id;
    const cs_real_t  *_b = gkb->b_tilda + s_id;
    const cs_lnum_t  block_size = CS_SBLOCK_BLOCK_SIZE;
    const cs_lnum_t  n_blocks = (n + block_size - 1) / block_size;
    const cs_lnum_t  n_sblocks = (n_blocks > 3) ? sqrt(n_blocks) : 1;
    const cs_lnum_t  blocks_in_sblocks =
      (n + block_size*n_sblocks - 1) / (block_size*n_sblocks);

    cs_real_t  *_q = gkb->q + s_id;
    cs_lnum_t  shift = 0;

    for (cs_lnum_t s = 0; s < n_sblocks; s++) { /* Loop on slices */

      double  s_beta2 = 0.0;

      for (cs_lnum_t b_id = 0; b_id < blocks_in_sblocks; b_id++) {

        const cs_lnum_t  start_id = shift;
        shift += block_size;
        if (shift > n)
          shift = n, b_id = blocks_in_sblocks;
        const cs_lnum_t  end_id = shift;

        double  _beta2 = 0.0;
        for (cs_lnum_t j = start_id; j < end_id; j++) {

          const  cs_real_t  b_ov_w = _b[j]/_w[j];
          _beta2 += b_ov_w*_b[j];
          _q[j] = b_ov_w;

        } /* Loop on block_size */

        s_beta2 += _beta2;

      } /* Loop on blocks */

      beta2 += s_beta2;

    } /* Loop on super-blocks */

  } /* OpenMP block */

  /* Parallel synchronization */

  cs_parall_sum(1, CS_DOUBLE, &beta2);

  /* Keep the value of beta = ||b||_{N^-1} */

  assert(beta2 > -DBL_MIN);
  gkb->beta = sqrt(beta2);

  /* Store M^-1.(b_f + gamma. Bt.N^-1.b_c) in b_tilda which is not useful
   * anymore */

  memcpy(gkb->b_tilda, gkb->v, gkb->n_u_dofs*sizeof(cs_real_t));

  if (fabs(gkb->beta) > FLT_MIN) {
    const  cs_real_t  inv_beta = 1./gkb->beta;
# pragma omp parallel for if (size > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < size; i++)
      gkb->q[i] *= inv_beta;
  }
  else {
    gkb->algo->cvg = CS_SLES_CONVERGED;
    return;
  }

  /* Solve M.w = Dt.q */

  _apply_div_op_transpose(div_op, gkb->q, gkb->dt_q);

  if (rset->ifs != NULL)
    cs_interface_set_sum(rset->ifs,
                         /* n_elts, stride, interlaced */
                         gkb->n_u_dofs, 1, false, CS_REAL_TYPE,
                         gkb->dt_q);

  cs_real_t  normalization = 1.0; /* TODO */

  gkb->algo->n_inner_iter
    += (gkb->algo->last_inner_iter =
        cs_cdo_solve_scalar_system(gkb->n_u_dofs,
                                   slesp,
                                   matrix,
                                   rset,
                                   normalization,
                                   false, /* rhs_redux */
                                   sles,
                                   gkb->v,
                                   gkb->dt_q));

  gkb->alpha = _face_gdot(rset, gkb->n_u_dofs, gkb->v, gkb->dt_q);
  assert(gkb->alpha > -DBL_MIN);
  gkb->alpha = sqrt(gkb->alpha);

  const double ov_alpha = 1./gkb->alpha;

  gkb->zeta = gkb->beta * ov_alpha;

  /* Initialize auxiliary vectors and first update of the solution vectors */

# pragma omp parallel for if (gkb->n_u_dofs > CS_THR_MIN)
  for (cs_lnum_t iu = 0; iu < gkb->n_u_dofs; iu++) {
    gkb->v[iu] *= ov_alpha;
    gkb->u_tilda[iu] = gkb->zeta * gkb->v[iu];
    gkb->m__v[iu] = ov_alpha * gkb->dt_q[iu];
  }

# pragma omp parallel for if (gkb->n_p_dofs > CS_THR_MIN)
  for (cs_lnum_t ip = 0; ip < gkb->n_p_dofs; ip++) {
    gkb->d[ip] = gkb->q[ip] * ov_alpha;
    p_c[ip] = -gkb->zeta * gkb->d[ip];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Test if one needs one more GKB iteration
 *
 * \param[in, out] gkb     pointer to a GKB builder structure
 */
/*----------------------------------------------------------------------------*/

static void
_gkb_cvg_test(cs_gkb_builder_t           *gkb)
{
  /* Update the sum of square of zeta values (used for renormalization) */

  cs_real_t  z2 = gkb->zeta*gkb->zeta;

  gkb->zeta_square_sum += z2;
  gkb->zeta_array[gkb->algo->n_algo_iter % gkb->z_size] = z2;

  /* Increment the number of Picard iterations */

  gkb->algo->n_algo_iter += 1;

  /* Compute the relative energy norm. The normalization arises from an
     iterative estimation of the initial error in the energy norm */

  const cs_real_t  prev_res = gkb->algo->res;

  int  n = gkb->z_size;
  if (gkb->algo->n_algo_iter < gkb->z_size)
    n = gkb->algo->n_algo_iter;

  cs_real_t  err2_energy = 0.;
  for (int i = 0; i < n; i++)
    err2_energy += gkb->zeta_array[i];

  double  tau = (gkb->gamma > 0) ?
    gkb->gamma*gkb->algo->param.rtol : gkb->algo->param.rtol;

  gkb->algo->res = sqrt(err2_energy);

  if (gkb->algo->n_algo_iter < 2)
    gkb->algo->res0 = gkb->algo->res;

  /* Set the convergence status */

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_SLES_DBG > 0
  cs_log_printf(CS_LOG_DEFAULT,
                "\nGKB.It%02d-- err2 = %6.4e ?<? tau * square_sum %6.4e\n",
                gkb->algo->n_algo_iter, err2_energy, tau*gkb->zeta_square_sum);
#endif

  if (err2_energy < tau * gkb->zeta_square_sum)
    gkb->algo->cvg = CS_SLES_CONVERGED;

  else if (gkb->algo->n_algo_iter >= gkb->algo->param.n_max_algo_iter)
    gkb->algo->cvg = CS_SLES_MAX_ITERATION;

  else if (gkb->algo->res > gkb->algo->param.dtol * prev_res ||
           gkb->algo->res > gkb->algo->param.dtol * gkb->algo->res0)
    gkb->algo->cvg = CS_SLES_DIVERGED;

  else
    gkb->algo->cvg = CS_SLES_ITERATING;

  if (gkb->algo->param.verbosity > 0)
    cs_log_printf(CS_LOG_DEFAULT,
                  "### GKB.It%d %5.3e %5d %6d z2:%6.4e renorm:%6.4e cvg:%d\n",
                  gkb->algo->n_algo_iter, gkb->algo->res,
                  gkb->algo->last_inner_iter, gkb->algo->n_inner_iter,
                  z2, sqrt(gkb->zeta_square_sum), gkb->algo->cvg);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Test if one needs one more Uzawa iteration in case of an Uzawa
 *         CG (conjugate gradient variant). The residual criterion has to be
 *         computed before calling this function.
 *
 * \param[in, out] uza     pointer to a Uzawa builder structure
 *
 * \return true (one more iteration) otherwise false
 */
/*----------------------------------------------------------------------------*/

static bool
_uza_cg_cvg_test(cs_uza_builder_t           *uza)
{
  /* Increment the number of algo. iterations */

  uza->algo->n_algo_iter += 1;

  /* Compute the new residual based on the norm of the divergence constraint */

  const cs_real_t  prev_res = uza->algo->res;
  const double  tau = fmax(uza->algo->param.rtol * uza->algo->normalization,
                           uza->algo->param.atol);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_SLES_DBG > 0
  cs_log_printf(CS_LOG_DEFAULT,
                "\nUZA-CG.It%02d-- res = %6.4e ?<? eps %6.4e\n",
                uza->algo->n_algo_iter, uza->algo->res, uza->algo->param.rtol);
#endif

  /* Set the convergence status */

  if (uza->algo->res < tau)
    uza->algo->cvg = CS_SLES_CONVERGED;

  else if (uza->algo->n_algo_iter >= uza->algo->param.n_max_algo_iter)
    uza->algo->cvg = CS_SLES_MAX_ITERATION;

  else if (uza->algo->res > uza->algo->param.dtol * prev_res ||
           uza->algo->res > uza->algo->param.dtol * uza->algo->res0)
    uza->algo->cvg = CS_SLES_DIVERGED;

  else
    uza->algo->cvg = CS_SLES_ITERATING;

  if (uza->algo->param.verbosity > 0)
    cs_log_printf(CS_LOG_DEFAULT,
                  "<UZACG.It%02d> res %5.3e | %4d %6d cvg%d | fit.eps %5.3e\n",
                  uza->algo->n_algo_iter, uza->algo->res,
                  uza->algo->last_inner_iter, uza->algo->n_inner_iter,
                  uza->algo->cvg, tau);

  if (uza->algo->cvg == CS_SLES_ITERATING)
    return true;
  else
    return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Test if one needs one more Uzawa iteration
 *
 * \param[in, out] uza     pointer to a Uzawa builder structure
 */
/*----------------------------------------------------------------------------*/

static void
_uza_cvg_test(cs_uza_builder_t           *uza)
{
  /* Increment the number of algo. iterations */

  uza->algo->n_algo_iter += 1;

  /* Compute the new residual based on the norm of the divergence constraint */

  const cs_real_t  prev_res = uza->algo->res;

  cs_real_t  res_square = cs_dot_wxx(uza->n_p_dofs, uza->inv_mp, uza->res_p);
  cs_parall_sum(1, CS_DOUBLE, &res_square);
  assert(res_square > -DBL_MIN);
  uza->algo->res = sqrt(res_square);

  double  tau = (uza->gamma > 0) ?
    uza->algo->param.rtol/sqrt(uza->gamma) : uza->algo->param.rtol;

  /* Set the convergence status */

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_SLES_DBG > 0
  cs_log_printf(CS_LOG_DEFAULT, "\nUZA.It%02d-- res = %6.4e ?<? eps %6.4e\n",
                uza->algo->n_algo_iter, uza->algo->res, uza->algo->param.rtol);
#endif

  if (uza->algo->res < tau)
    uza->algo->cvg = CS_SLES_CONVERGED;

  else if (uza->algo->n_algo_iter >= uza->algo->param.n_max_algo_iter)
    uza->algo->cvg = CS_SLES_MAX_ITERATION;

  else if (uza->algo->res > uza->algo->param.dtol * prev_res)
    uza->algo->cvg = CS_SLES_DIVERGED;

  else
    uza->algo->cvg = CS_SLES_ITERATING;

  if (uza->algo->param.verbosity > 0)
    cs_log_printf(CS_LOG_DEFAULT, "### UZA.It%02d-- %5.3e %5d %6d cvg:%d\n",
                  uza->algo->n_algo_iter, uza->algo->res,
                  uza->algo->last_inner_iter, uza->algo->n_inner_iter,
                  uza->algo->cvg);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Test if one needs one more Uzawa iteration in case of an incremental
 *         formulation
 *
 * \param[in]      delta_u_l2   value of the weighted L2 norm of delta_u
 * \param[in, out] uza          pointer to a Uzawa builder structure
 *
 * \return true if one more iteration is needed otherwise false
 */
/*----------------------------------------------------------------------------*/

static bool
_uza_incr_cvg_test(cs_real_t                   delta_u_l2,
                   cs_uza_builder_t           *uza)
{
  cs_iter_algo_t  *algo = uza->algo;

  /* Compute the new residual based on the norm of the divergence constraint */

  algo->prev_res = algo->res;

  cs_real_t  res_square = cs_dot_wxx(uza->n_p_dofs, uza->inv_mp, uza->d__v);
  cs_parall_sum(1, CS_DOUBLE, &res_square);
  assert(res_square > -DBL_MIN);
  cs_real_t  divu_l2 = sqrt(res_square);

  algo->res = fmax(delta_u_l2, divu_l2);

  if (algo->n_algo_iter == 0) { /* First call */

    algo->res0 = algo->res;
    algo->normalization = algo->res0;

    if (algo->param.verbosity > 1)
      cs_log_printf(CS_LOG_DEFAULT, "### UZAi.res0:%5.3e modified rtol:%5.3e"
                    " atol:%5.3e\n",
                    algo->res0, algo->normalization*algo->param.rtol,
                    algo->param.atol);

  }

  /* Update the convergence status */

  cs_iter_algo_update_cvg(algo);

  if (algo->param.verbosity > 0)
    cs_log_printf(CS_LOG_DEFAULT,
                  "### UZAi.It%02d %5.3e %5d %6d cvg:%d div:%5.3e, du:%5.3e\n",
                  algo->n_algo_iter, algo->res,
                  algo->last_inner_iter, algo->n_inner_iter,
                  algo->cvg, divu_l2, delta_u_l2);

  if (algo->cvg == CS_SLES_ITERATING)
    return true;
  else
    return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check and possibly fix the settings related to the Schur complement
 *         in case of AMG preconditioner
 *
 * \param[in, out] schur_slesp   pointer to the SLES parameters for the Schur
 */
/*----------------------------------------------------------------------------*/

static void
_set_schur_sles(cs_param_sles_t   *schur_slesp)
{
  if (schur_slesp->precond == CS_PARAM_PRECOND_AMG) {

    if (schur_slesp->amg_type == CS_PARAM_AMG_NONE) {

      cs_param_sles_class_t  ret_class =
        cs_param_sles_check_class(CS_PARAM_SLES_CLASS_HYPRE);

      if (ret_class == CS_PARAM_SLES_N_CLASSES)
        schur_slesp->amg_type = CS_PARAM_AMG_HOUSE_K;
      else
        schur_slesp->amg_type = CS_PARAM_AMG_HYPRE_BOOMER_V;

    }
    else {

      cs_param_sles_class_t  ret_class =
        cs_param_sles_get_class_from_amg(schur_slesp->amg_type);

      /* Modify the default settings if needed */

      if (ret_class != schur_slesp->solver_class &&
          schur_slesp->solver_class == CS_PARAM_SLES_CLASS_CS)
        schur_slesp->solver_class = ret_class;

    }

    cs_param_sles_check_amg(schur_slesp);

  } /* Check AMG settings */

  int ier = cs_param_sles_set(false, schur_slesp);

  if (ier == -1)
    bft_error(__FILE__, __LINE__, 0,
              "%s: The requested class of solvers is not available"
              " for the system %s\n Please modify your settings.",
              __func__, schur_slesp->name);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create an empty cs_cdofb_monolithic_sles_t structure
 *
 * \param[in] n_faces     number of faces (interior + border)
 * \param[in] n_cells     number of cells
 *
 * \return a pointer to a newly allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_cdofb_monolithic_sles_t *
cs_cdofb_monolithic_sles_create(cs_lnum_t    n_faces,
                                cs_lnum_t    n_cells)
{
  cs_cdofb_monolithic_sles_t  *msles = NULL;

  BFT_MALLOC(msles, 1, cs_cdofb_monolithic_sles_t);

  msles->div_op = NULL;

  msles->graddiv_coef = 0.;

  msles->sles = NULL;
  msles->schur_sles = NULL;

  msles->n_faces = n_faces;
  msles->n_cells = n_cells;

  msles->u_f = NULL;
  msles->p_c = NULL;

  return msles;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a part of the structure
 *
 * \param[in, out]  msles   pointer to the structure to clean
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_monolithic_sles_clean(cs_cdofb_monolithic_sles_t   *msles)
{
  if (msles == NULL)
    return;

  cs_sles_free(msles->sles);
  cs_sles_free(msles->schur_sles);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free memory related to cs_cdofb_monolithic_sles_t structure
 *
 * \param[in, out]  p_msles  double pointer to the structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_monolithic_sles_free(cs_cdofb_monolithic_sles_t   **p_msles)
{
  cs_cdofb_monolithic_sles_t  *msles = *p_msles;

  if (msles == NULL)
    return;

  /* sles are freed elsewhere */

  BFT_FREE(msles->div_op);

  /* Other pointer are shared, thus no free at this stage */

  BFT_FREE(msles);
  *p_msles = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set pointers to shared structures
 *
 * \param[in]  connect  pointer to cdo connectivities
 * \param[in]  quant    pointer to additional mesh quantities
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_monolithic_sles_init_sharing(const cs_cdo_connect_t        *connect,
                                      const cs_cdo_quantities_t     *quant)
{
  /* Assign static const pointers */

  cs_shared_connect = connect;
  cs_shared_quant = quant;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free if needed structure(s) associated CDO face-based schemes with
 *         a monolithic velocity-pressure coupling
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_monolithic_sles_finalize(void)
{
  if (_petsc_hook_context != NULL)
    BFT_FREE(_petsc_hook_context); /* contains only shared pointers */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Start setting-up the Navier-Stokes equations when a monolithic
 *         algorithm is used to couple the system.
 *         No mesh information is available at this stage.
 *         nsp is not declared as const to avoid compilation warnings but
 *         it should be modified at this stage.
 *
 * \param[in]      nsp      pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_monolithic_set_sles(cs_navsto_param_t    *nsp,
                             void                 *context)
{
  cs_cdofb_monolithic_t  *sc = context;
  cs_navsto_monolithic_t *cc = (cs_navsto_monolithic_t *)sc->coupling_context;
  cs_equation_t  *mom_eq = cc->momentum;

  const cs_navsto_param_sles_t  *nslesp = nsp->sles_param;
  cs_equation_param_t  *mom_eqp = cs_equation_get_param(mom_eq);
  cs_param_sles_t  *mom_slesp = mom_eqp->sles_param;
  int  field_id = cs_equation_get_field_id(mom_eq);

  mom_slesp->field_id = field_id;
  if (mom_slesp->amg_type == CS_PARAM_AMG_NONE) {

    cs_param_sles_class_t  ret_class =
      cs_param_sles_check_class(CS_PARAM_SLES_CLASS_HYPRE);

    if (ret_class == CS_PARAM_SLES_N_CLASSES)
      mom_slesp->amg_type = CS_PARAM_AMG_HOUSE_K;
    else
      mom_slesp->amg_type = CS_PARAM_AMG_HYPRE_BOOMER_V;

  }

  switch (nslesp->strategy) {

  case CS_NAVSTO_SLES_EQ_WITHOUT_BLOCK: /* "Classical" way to set SLES */
    cs_equation_param_set_sles(mom_eqp);
    break;

  case CS_NAVSTO_SLES_GKB_SATURNE:
    /* Set solver and preconditioner for solving M = A + zeta * Bt*N^-1*B
     * Notice that zeta can be equal to 0 */
    cs_equation_param_set_sles(mom_eqp);
    break;

  case CS_NAVSTO_SLES_MINRES:
  case CS_NAVSTO_SLES_GCR:
    cs_equation_param_set_sles(mom_eqp);
    break;

  case CS_NAVSTO_SLES_DIAG_SCHUR_MINRES:
  case CS_NAVSTO_SLES_DIAG_SCHUR_GCR:
  case CS_NAVSTO_SLES_LOWER_SCHUR_GCR:
  case CS_NAVSTO_SLES_SGS_SCHUR_GCR:
  case CS_NAVSTO_SLES_UPPER_SCHUR_GCR:
  case CS_NAVSTO_SLES_UZAWA_SCHUR_GCR:
    {
      /* Set solver and preconditioner for solving A */

      cs_equation_param_set_sles(mom_eqp);

      /* Set the solver for the compatible Laplacian (the related SLES is
         defined using the system name instead of the field id since this is an
         auxiliary system) */

      _set_schur_sles(nslesp->schur_sles_param);
    }
    break;

  case CS_NAVSTO_SLES_UZAWA_CG:
    {
      /* Set solver and preconditioner for solving A */

      cs_equation_param_set_sles(mom_eqp);

      /* Set the solver for the compatible Laplacian (the related SLES is
         defined using the system name instead of the field id since this is an
         auxiliary system) */

      _set_schur_sles(nslesp->schur_sles_param);
    }
    break;

  case CS_NAVSTO_SLES_UZAWA_AL:
     /* Set solver and preconditioner for solving M = A + zeta * Bt*N^-1*B
      * Notice that zeta can be equal to 0 */

    cs_equation_param_set_sles(mom_eqp);
    break;

#if defined(HAVE_PETSC)
  /* Strategies available before the 3.11 version of PETSc */

  case CS_NAVSTO_SLES_ADDITIVE_GMRES_BY_BLOCK:
    cs_sles_petsc_init();
    _initialize_petsc_hook_context(nsp, sc);
    cs_sles_petsc_define(field_id,
                         NULL,
                         MATMPIAIJ,
                         _additive_amg_hook,
                         (void *)_petsc_hook_context);
    break;

  case CS_NAVSTO_SLES_DIAG_SCHUR_GMRES:
    cs_sles_petsc_init();
    _initialize_petsc_hook_context(nsp, sc);
    cs_sles_petsc_define(field_id,
                         NULL,
                         MATMPIAIJ,
                         _diag_schur_hook,
                         (void *)_petsc_hook_context);
    break;

  case CS_NAVSTO_SLES_MULTIPLICATIVE_GMRES_BY_BLOCK:
    cs_sles_petsc_init();
    _initialize_petsc_hook_context(nsp, sc);
    cs_sles_petsc_define(field_id,
                         NULL,
                         MATMPIAIJ,
                         _multiplicative_hook,
                         (void *)_petsc_hook_context);
    break;

  case CS_NAVSTO_SLES_NOTAY_TRANSFORM:
    cs_sles_petsc_init();
    _initialize_petsc_hook_context(nsp, sc);
    cs_sles_petsc_define(field_id,
                         NULL,
                         MATMPIAIJ,
                         _notay_hook,
                         (void *)_petsc_hook_context);
    break;

  case CS_NAVSTO_SLES_UPPER_SCHUR_GMRES:
    cs_sles_petsc_init();
    _initialize_petsc_hook_context(nsp, sc);
    cs_sles_petsc_define(field_id,
                         NULL,
                         MATMPIAIJ,
                         _upper_schur_hook,
                         (void *)_petsc_hook_context);
    break;

  /* Golub-Kahan Bi-diagonalization is available starting from the 3.11 version
     of PETSc */
#if PETSC_VERSION_GE(3,11,0)
  case CS_NAVSTO_SLES_GKB_PETSC:
    cs_sles_petsc_init();
    _initialize_petsc_hook_context(nsp, sc);
    cs_sles_petsc_define(field_id,
                         NULL,
                         MATMPIAIJ,
                         _gkb_hook,
                         (void *)_petsc_hook_context);
    break;

  case CS_NAVSTO_SLES_GKB_GMRES:
    cs_sles_petsc_init();
    _initialize_petsc_hook_context(nsp, sc);
    cs_sles_petsc_define(field_id,
                         NULL,
                         MATMPIAIJ,
                         _gkb_precond_hook,
                         (void *)_petsc_hook_context);
    break;
#else  /* PETSC version < 3.11 */
  case CS_NAVSTO_SLES_GKB_PETSC:
  case CS_NAVSTO_SLES_GKB_GMRES:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid strategy for solving the linear system %s\n"
              " PETSc 3.11.x or greater is required with this option.\n",
              __func__, mom_eqp->name);
    break;
#endif /* PETSc version */

#else  /* no HAVE_PETSC */

  case CS_NAVSTO_SLES_ADDITIVE_GMRES_BY_BLOCK:
  case CS_NAVSTO_SLES_DIAG_SCHUR_GMRES:
  case CS_NAVSTO_SLES_GKB_PETSC:
  case CS_NAVSTO_SLES_GKB_GMRES:
  case CS_NAVSTO_SLES_MULTIPLICATIVE_GMRES_BY_BLOCK:
  case CS_NAVSTO_SLES_UPPER_SCHUR_GMRES:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid strategy for solving the linear system %s\n"
              " PETSc is required with this option.\n"
              " Please use a version of code_saturne built with PETSc.",
              __func__, mom_eqp->name);
    break;

#endif  /* HAVE_PETSC */

  case CS_NAVSTO_SLES_MUMPS:
#if defined(HAVE_MUMPS)
    if (mom_slesp->solver != CS_PARAM_ITSOL_MUMPS &&
        mom_slesp->solver != CS_PARAM_ITSOL_MUMPS_LDLT &&
        mom_slesp->solver != CS_PARAM_ITSOL_MUMPS_FLOAT &&
        mom_slesp->solver != CS_PARAM_ITSOL_MUMPS_FLOAT_LDLT)
      mom_slesp->solver = CS_PARAM_ITSOL_MUMPS;

    cs_sles_mumps_define(field_id,
                         NULL,
                         mom_slesp,
                         cs_user_sles_mumps_hook,
                         NULL);
#else
#if defined(HAVE_PETSC)
#if defined(PETSC_HAVE_MUMPS)
    cs_sles_petsc_init();
    cs_sles_petsc_define(field_id,
                         NULL,
                         MATMPIAIJ,
                         _mumps_hook,
                         (void *)mom_eqp);
#else
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid strategy for solving the linear system %s\n"
              " PETSc with MUMPS is required with this option.\n",
              __func__, mom_eqp->name);
#endif  /* PETSC_HAVE_MUMPS */
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid strategy for solving the linear system %s\n"
              " Neither PETSc nor MUMPS is available.\n",
              __func__, mom_eqp->name);
#endif  /* HAVE_PETSC */
#endif  /* HAVE_MUMPS */
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid strategy for solving the linear system %s\n",
              __func__, mom_eqp->name);
  }

  /* Define the level of verbosity for SLES structure */

  if (mom_slesp->verbosity > 1) {

    cs_sles_t  *sles = cs_sles_find_or_add(field_id, NULL);

    /* Set verbosity */

    cs_sles_set_verbosity(sles, mom_slesp->verbosity);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve a linear system arising from the discretization of the
 *         Navier-Stokes equation with a CDO face-based approach.
 *         The full system is treated as one block and solved as it is.
 *         In this situation, PETSc or MUMPS are usually considered.
 *
 * \param[in]      nsp      pointer to a cs_navsto_param_t structure
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in]      sh       pointer to a cs_cdo_system_helper_t structure
 * \param[in, out] slesp    pointer to a set of parameters to drive the SLES
 * \param[in, out] msles    pointer to a cs_cdofb_monolithic_sles_t structure
 *
 * \return the (cumulated) number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdofb_monolithic_solve(const cs_navsto_param_t       *nsp,
                          const cs_equation_param_t     *eqp,
                          const cs_cdo_system_helper_t  *sh,
                          cs_param_sles_t               *slesp,
                          cs_cdofb_monolithic_sles_t    *msles)
{
  assert(sh != NULL);
  assert(sh->n_blocks == 1);

  const cs_range_set_t  *range_set = cs_cdo_system_get_range_set(sh, 0);
  const cs_matrix_t  *matrix = cs_cdo_system_get_matrix(sh, 0);
  const cs_lnum_t  n_cols = cs_matrix_get_n_columns(matrix);
  const cs_lnum_t  n_faces = msles->n_faces;
  const cs_lnum_t  n_cells = msles->n_cells;
  const cs_lnum_t  n_scatter_elts = 3*n_faces + n_cells;
  const cs_navsto_param_sles_t  *nslesp = nsp->sles_param;

  /* De-interlace the velocity array and the rhs for the face DoFs */

  cs_real_t  *sol = NULL;
  BFT_MALLOC(sol, CS_MAX(n_cols, n_scatter_elts), cs_real_t);

  cs_real_t  *b = NULL;
  BFT_MALLOC(b, n_scatter_elts, cs_real_t);

# pragma omp parallel for if (CS_THR_MIN > n_faces)     \
  shared(msles, sol, b) firstprivate(n_faces)
  for (cs_lnum_t f = 0; f < n_faces; f++) {

    sol[f            ] = msles->u_f[3*f];
    sol[f +   n_faces] = msles->u_f[3*f+1];
    sol[f + 2*n_faces] = msles->u_f[3*f+2];

    b[f            ] = sh->rhs[3*f];
    b[f +   n_faces] = sh->rhs[3*f+1];
    b[f + 2*n_faces] = sh->rhs[3*f+2];

  }

  /* Add the pressure related elements */

  memcpy(sol + 3*n_faces, msles->p_c, n_cells*sizeof(cs_real_t));
  memcpy(b + 3*n_faces, sh->rhs + 3*n_faces, n_cells*sizeof(cs_real_t));

  if (nslesp->strategy == CS_NAVSTO_SLES_NOTAY_TRANSFORM) {

# pragma omp parallel for if (CS_THR_MIN > n_cells)     \
  shared(b) firstprivate(n_faces)
    for (cs_lnum_t f = 3*n_faces; f < n_scatter_elts; f++)
      b[f] = -1.0*b[f];

  }

  int  n_iters = 0;
  double  residual = DBL_MAX;

  /* Prepare solving (handle parallelism) */

  cs_cdo_solve_prepare_system(1,     /* stride */
                              false, /* interlace (managed here) */
                              n_scatter_elts,
                              range_set,
                              true,  /* rhs_redux */
                              sol, b);

  /* Solve the linear solver */

  const double  r_norm = 1.0; /* No renormalization by default (TODO) */

  cs_real_t  rtol = slesp->eps;

  if (nslesp->strategy == CS_NAVSTO_SLES_UPPER_SCHUR_GMRES              ||
      nslesp->strategy == CS_NAVSTO_SLES_DIAG_SCHUR_GMRES               ||
      nslesp->strategy == CS_NAVSTO_SLES_MULTIPLICATIVE_GMRES_BY_BLOCK  ||
      nslesp->strategy == CS_NAVSTO_SLES_NOTAY_TRANSFORM                ||
      nslesp->strategy == CS_NAVSTO_SLES_ADDITIVE_GMRES_BY_BLOCK)
    rtol = nslesp->il_algo_param.rtol;

  cs_sles_convergence_state_t  code = cs_sles_solve(msles->sles,
                                                    matrix,
                                                    rtol,
                                                    r_norm,
                                                    &n_iters,
                                                    &residual,
                                                    b,
                                                    sol,
                                                    0,      /* aux. size */
                                                    NULL);  /* aux. buffers */

  /* Output information about the convergence of the resolution */

  if (slesp->verbosity > 1)
    cs_log_printf(CS_LOG_DEFAULT, "  <%20s/sles_cvg_code=%-d> n_iters %d |"
                  " residual % -8.4e\n",
                  eqp->name, code, n_iters, residual);

  /* sol is computed and stored in a "gather" view. Switch to a "scatter"
     view */

  cs_range_set_scatter(range_set,
                       CS_REAL_TYPE, 1, /* type and stride */
                       sol, sol);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_SLES_DBG > 1
  cs_range_set_scatter(range_set,
                       CS_REAL_TYPE, 1, /* type and stride */
                       b, b);

  cs_dbg_fprintf_system(eqp->name,
                        -1,
                        CS_CDOFB_MONOLITHIC_SLES_DBG,
                        sol, b, 3*n_faces);
#endif

  /* Switch from sol (not interlaced) to u_f and p_c */

  cs_real_t  *u_f = msles->u_f;

  /* Copy the part of the solution array related to the pressure in cells */

  memcpy(msles->p_c, sol + 3*n_faces, n_cells*sizeof(cs_real_t));

  if (nslesp->strategy == CS_NAVSTO_SLES_NOTAY_TRANSFORM) {

    cs_real_t  *grad_p = NULL, *mat_diag = NULL;

    /* Compute the pressure gradient */

    BFT_MALLOC(grad_p, 3*n_faces, cs_real_t);

    _apply_div_op_transpose(msles->div_op, msles->p_c, grad_p);

    /* grad_p is build cellwise. Perform the parallel synchronization. */

    if (range_set->ifs != NULL)
      cs_interface_set_sum(range_set->ifs,
                         /* n_elts, stride, interlaced */
                           3*n_faces, 1, false, CS_REAL_TYPE,
                           grad_p);

    /* Retrieve the diagonal of the matrix in a "scatter" view */

    BFT_MALLOC(mat_diag, n_scatter_elts, cs_real_t);

    /* diag is stored in a "gather view". Switch to a "scatter view" to make
       the change of variable */

    cs_range_set_scatter(range_set,
                         CS_REAL_TYPE,
                         1,         /* treated as scalar-valued up to now */
                         cs_matrix_get_diagonal(matrix), /* gathered view */
                         mat_diag);                      /* scatter view */

    const double  alpha = cs_navsto_param_get_notay_scaling();
    const cs_real_t  *dx = mat_diag, *dy = mat_diag + n_faces;
    const cs_real_t  *dz = mat_diag + 2*n_faces;
    const cs_real_t  *solx = sol, *soly = sol+n_faces, *solz = sol+2*n_faces;

# pragma omp parallel for if (CS_THR_MIN > n_faces)                     \
  shared(dx, dy, dz, solx, soly, solz) firstprivate(n_faces)
    for (cs_lnum_t f = 0; f < n_faces; f++) {
      u_f[3*f  ] = solx[f] - alpha * grad_p[3*f]/dx[f];
      u_f[3*f+1] = soly[f] - alpha * grad_p[3*f+1]/dy[f];
      u_f[3*f+2] = solz[f] - alpha * grad_p[3*f+2]/dz[f];
    }

    BFT_FREE(grad_p);
    BFT_FREE(mat_diag);

  }
  else { /* Other strategies: No change of variable */

# pragma omp parallel for if (CS_THR_MIN > n_faces)     \
  shared(msles, sol) firstprivate(n_faces)
    for (cs_lnum_t f = 0; f < n_faces; f++) {
      u_f[3*f  ] = sol[f          ];
      u_f[3*f+1] = sol[f+  n_faces];
      u_f[3*f+2] = sol[f+2*n_faces];
    }

  }

  /* Free what can be freed at this stage */

  BFT_FREE(sol);
  BFT_FREE(b);

  return n_iters;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve a linear system arising from the discretization of the
 *        Navier-Stokes equation with a CDO face-based approach. The system is
 *        split into blocks to enable more efficient preconditioning
 *        techniques. The main iterative solver is a Krylov solver such as GCR,
 *        or MINRES
 *
 * \param[in]      nsp      pointer to a cs_navsto_param_t structure
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in]      sh       pointer to a cs_cdo_system_helper_t structure
 * \param[in, out] slesp    pointer to a set of parameters to drive the SLES
 * \param[in, out] msles    pointer to a cs_cdofb_monolithic_sles_t structure
 *
 * \return the (cumulated) number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdofb_monolithic_krylov_block_precond(const cs_navsto_param_t       *nsp,
                                         const cs_equation_param_t     *eqp,
                                         const cs_cdo_system_helper_t  *sh,
                                         cs_param_sles_t               *slesp,
                                         cs_cdofb_monolithic_sles_t    *msles)
{
  CS_UNUSED(eqp);

  if (msles == NULL)
    return 0;

  if (sh->type != CS_CDO_SYSTEM_SADDLE_POINT)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid type of system: saddle-point system expected\n",
              __func__);
  if (sh->n_col_blocks != 2)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid type of system: 2x2 saddle-point system expected\n",
              __func__);

  cs_cdo_system_block_t  *vel_block = sh->blocks[0];
  cs_cdo_system_block_t  *div_block = sh->blocks[1];

  if (vel_block->type != CS_CDO_SYSTEM_BLOCK_DEFAULT)
    bft_error(__FILE__, __LINE__, 0,
              "%s: A default block for the velocity is expected.",
              __func__);

  assert(sh->col_block_sizes[0] == 3*msles->n_faces);
  assert(sh->col_block_sizes[1] == msles->n_cells);
  assert(div_block->info.stride == 3);

  /* Set the structure to manage the iterative solver */

  const cs_navsto_param_sles_t  *nslesp = nsp->sles_param;

  cs_iter_algo_t  *saddle_algo = cs_iter_algo_create(nslesp->il_algo_param);

  /* Set the saddle-point system */
  /* --------------------------- */

  cs_cdo_system_dblock_t  *vel_dblock = vel_block->block_pointer;
  cs_saddle_system_t  *ssys = NULL;

  BFT_MALLOC(ssys, 1, cs_saddle_system_t);

  ssys->n_m11_matrices = 1;
  BFT_MALLOC(ssys->m11_matrices, 1, cs_matrix_t *);
  ssys->m11_matrices[0] = vel_dblock->matrix;
  ssys->x1_size = sh->col_block_sizes[0];
  ssys->max_x1_size =
    CS_MAX(ssys->x1_size, cs_matrix_get_n_columns(vel_dblock->matrix));
  ssys->rhs1 = sh->rhs_array[0];

  cs_cdo_system_ublock_t  *div_ublock = div_block->block_pointer;
  assert(msles->div_op == div_ublock->values);

  ssys->x2_size = sh->col_block_sizes[1];
  ssys->max_x2_size = cs_glob_mesh->n_cells_with_ghosts;
  ssys->rhs2 = sh->rhs_array[1];

  ssys->m21_stride = 3;
  ssys->m21_unassembled = msles->div_op;
  ssys->m21_adjacency = div_ublock->adjacency;

  ssys->rset = vel_dblock->range_set;

  /* u_f is allocated to 3*n_faces (the size of the scatter view but during the
     resolution process one need a vector at least of size n_cols of the matrix
     m11. */

  cs_real_t  *xu = NULL;
  BFT_MALLOC(xu, ssys->max_x1_size, cs_real_t);
  memcpy(xu, msles->u_f, ssys->x1_size*sizeof(cs_real_t));

  switch (nslesp->strategy) {

  case CS_NAVSTO_SLES_DIAG_SCHUR_GCR:
  case CS_NAVSTO_SLES_LOWER_SCHUR_GCR:
  case CS_NAVSTO_SLES_SGS_SCHUR_GCR:
  case CS_NAVSTO_SLES_UPPER_SCHUR_GCR:
  case CS_NAVSTO_SLES_UZAWA_SCHUR_GCR:
    {
      /* Default */

      cs_param_precond_block_t
        preblock_type = CS_PARAM_PRECOND_BLOCK_UPPER_TRIANGULAR;

      if (nslesp->strategy == CS_NAVSTO_SLES_DIAG_SCHUR_GCR)
        preblock_type = CS_PARAM_PRECOND_BLOCK_DIAG;
      else if (nslesp->strategy == CS_NAVSTO_SLES_LOWER_SCHUR_GCR)
        preblock_type = CS_PARAM_PRECOND_BLOCK_LOWER_TRIANGULAR;
      else if (nslesp->strategy == CS_NAVSTO_SLES_SGS_SCHUR_GCR)
        preblock_type = CS_PARAM_PRECOND_BLOCK_SYM_GAUSS_SEIDEL;
      else if (nslesp->strategy == CS_NAVSTO_SLES_UZAWA_SCHUR_GCR)
        preblock_type = CS_PARAM_PRECOND_BLOCK_UZAWA;

      /* Define block preconditionning */

      cs_saddle_block_precond_t  *sbp =
        cs_saddle_block_precond_create(preblock_type,
                                       nslesp->schur_approximation,
                                       slesp,
                                       msles->sles);

      /* Define an approximation of the Schur complement */

      _schur_approximation(nsp, ssys, msles->schur_sles, sbp);

      /* Call the inner linear algorithm */

      cs_saddle_gcr(nslesp->il_algo_restart, ssys, sbp,
                    xu, msles->p_c, saddle_algo);

      cs_saddle_block_precond_free(&sbp);
    }
    break;

  case CS_NAVSTO_SLES_DIAG_SCHUR_MINRES:
    {
      /* Define block preconditionning */

      cs_saddle_block_precond_t  *sbp =
        cs_saddle_block_precond_create(CS_PARAM_PRECOND_BLOCK_DIAG,
                                       nslesp->schur_approximation,
                                       slesp,
                                       msles->sles);

      /* Define an approximation of the Schur complement */

      _schur_approximation(nsp, ssys, msles->schur_sles, sbp);

      /* Call the inner linear algorithm */

      cs_saddle_minres(ssys, sbp, xu, msles->p_c, saddle_algo);

      cs_saddle_block_precond_free(&sbp);
    }
    break;

  case CS_NAVSTO_SLES_GCR:   /* No block preconditioning */
    cs_saddle_gcr(nslesp->il_algo_restart, ssys, NULL,
                  xu, msles->p_c, saddle_algo);
    break;

  case CS_NAVSTO_SLES_MINRES:   /* No block preconditioning */
    cs_saddle_minres(ssys, NULL, xu, msles->p_c, saddle_algo);
    break;

  case CS_NAVSTO_SLES_USER:     /* User-defined strategy */
    {
      /* Define block preconditionning by default */

      cs_saddle_block_precond_t  *sbp =
        cs_saddle_block_precond_create(CS_PARAM_PRECOND_BLOCK_DIAG,
                                       nslesp->schur_approximation,
                                       slesp,
                                       msles->sles);

      cs_user_navsto_sles_solve(nslesp, ssys, sbp, xu, msles->p_c, saddle_algo);

      cs_saddle_block_precond_free(&sbp);
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid strategy to solve the system.\n"
              "Please used a Krylov-based iterative solver.", __func__);
    break;
  }

  memcpy(msles->u_f, xu, ssys->x1_size*sizeof(cs_real_t));

  /* Free the saddle-point system */

  BFT_FREE(xu);
  BFT_FREE(ssys->m11_matrices);
  BFT_FREE(ssys); /* only shared pointers inside */

  int n_algo_iter = saddle_algo->n_algo_iter;

  if (nslesp->il_algo_param.verbosity > 1)
    cs_log_printf(CS_LOG_DEFAULT,
                  " -cvg- inner_algo: cumulated_iters: %d\n",
                  saddle_algo->n_inner_iter);


  BFT_FREE(saddle_algo);

  return n_algo_iter;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Use the GKB algorithm to solve the saddle-point problem arising
 *         from CDO-Fb schemes for Stokes and Navier-Stokes with a monolithic
 *         coupling
 *
 * \param[in]      nsp      pointer to a cs_navsto_param_t structure
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in]      sh       pointer to a cs_cdo_system_helper_t structure
 * \param[in, out] slesp    pointer to a set of parameters to drive the SLES
 * \param[in, out] msles    pointer to a cs_cdofb_monolithic_sles_t structure
 *
 * \return the cumulated number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdofb_monolithic_gkb_solve(const cs_navsto_param_t       *nsp,
                              const cs_equation_param_t     *eqp,
                              const cs_cdo_system_helper_t  *sh,
                              cs_param_sles_t               *slesp,
                              cs_cdofb_monolithic_sles_t    *msles)
{
  assert(nsp != NULL);
  const cs_navsto_param_sles_t  *nslesp = nsp->sles_param;

  assert(nslesp->strategy == CS_NAVSTO_SLES_GKB_SATURNE);

  const cs_real_t  *vol = cs_shared_quant->cell_vol;
  const cs_real_t  *div_op = msles->div_op;
  const cs_real_t  gamma = msles->graddiv_coef;
  const cs_range_set_t  *range_set = cs_cdo_system_get_range_set(sh, 0);
  const cs_matrix_t  *matrix = cs_cdo_system_get_matrix(sh, 0);

  cs_real_t  *u_f = msles->u_f;
  cs_real_t  *p_c = msles->p_c;
  cs_real_t  *b_f = sh->rhs;
  cs_real_t  *b_c = sh->rhs + 3*msles->n_faces;

  /* Allocate and initialize the GKB builder structure */

  cs_gkb_builder_t  *gkb = _init_gkb_builder(nsp,
                                             gamma,
                                             3*msles->n_faces,
                                             msles->n_cells);

  /* Transformation of the initial saddle-point system */

  _transform_gkb_system(matrix, range_set, eqp, nslesp, div_op,
                        slesp, gkb, msles->sles, u_f, b_f, b_c);

  /* Initialization */

  _init_gkb_algo(matrix, range_set, div_op, slesp, gkb, msles->sles, p_c);

  /* Main loop */
  /* ========= */

  while (gkb->algo->cvg == CS_SLES_ITERATING) {

    /* Compute g (store as an update of d__v), q */

    _apply_div_op(div_op, gkb->v, gkb->d__v);

#   pragma omp parallel for if (gkb->n_p_dofs > CS_THR_MIN)
    for (cs_lnum_t ip = 0; ip < gkb->n_p_dofs; ip++) {
      gkb->d__v[ip] /= vol[ip];
      gkb->d__v[ip] -= gkb->alpha * gkb->q[ip];
    }

    /* Compute beta */

    gkb->beta = cs_dot_wxx(gkb->n_p_dofs, vol, gkb->d__v);
    cs_parall_sum(1, CS_DOUBLE, &(gkb->beta));
    assert(gkb->beta > -DBL_MIN);
    gkb->beta = sqrt(gkb->beta);

    const double  ov_beta = 1./gkb->beta;

#   pragma omp parallel for if (gkb->n_p_dofs > CS_THR_MIN)
    for (cs_lnum_t ip = 0; ip < gkb->n_p_dofs; ip++)
      gkb->q[ip] = ov_beta*gkb->d__v[ip];

    /* Solve M.w_tilda = Dt.q */

    _apply_div_op_transpose(div_op, gkb->q, gkb->dt_q);

    if (range_set->ifs != NULL)
      cs_interface_set_sum(range_set->ifs,
                           gkb->n_u_dofs,
                           1, false, CS_REAL_TYPE, /* stride, interlaced */
                           gkb->dt_q);

    /* Prepare update of m__v:
     *  m__v(k+1) = 1/alpha(k+1) * (dt_q - beta*m__v(k)) */

#   pragma omp parallel for if (gkb->n_u_dofs > CS_THR_MIN)
    for (cs_lnum_t iu = 0; iu < gkb->n_u_dofs; iu++) {
      gkb->m__v[iu] *= -gkb->beta;
      gkb->m__v[iu] +=  gkb->dt_q[iu];
    }

    cs_real_t  normalization = gkb->alpha; /* TODO */
    gkb->algo->n_inner_iter
      += (gkb->algo->last_inner_iter =
          cs_cdo_solve_scalar_system(gkb->n_u_dofs,
                                     slesp,
                                     matrix,
                                     range_set,
                                     normalization,
                                     false, /* rhs_redux */
                                     msles->sles,
                                     gkb->v,
                                     gkb->m__v));

    /* Compute alpha */

    gkb->alpha = _face_gdot(range_set, gkb->n_u_dofs, gkb->v, gkb->m__v);
    assert(gkb->alpha > -DBL_MIN);
    gkb->alpha = sqrt(gkb->alpha);

    const double ov_alpha = 1./gkb->alpha;

    /* zeta(k+1) = -beta/alpha * zeta(k) */

    gkb->zeta *= -gkb->beta * ov_alpha;

    /* Update vectors and solutions */

#   pragma omp parallel for if (gkb->n_u_dofs > CS_THR_MIN)
    for (cs_lnum_t iu = 0; iu < gkb->n_u_dofs; iu++) {
      gkb->v[iu] *= ov_alpha;
      gkb->u_tilda[iu] += gkb->zeta * gkb->v[iu];

      /* Last step: m__v(k+1) = 1/alpha(k+1) * (dt_q - beta*m__v(k)) */

      gkb->m__v[iu] *= ov_alpha;
    }

#   pragma omp parallel for if (gkb->n_p_dofs > CS_THR_MIN)
    for (cs_lnum_t ip = 0; ip < gkb->n_p_dofs; ip++) {
      gkb->d[ip] = ov_alpha * (gkb->q[ip] - gkb->beta*gkb->d[ip]);
      p_c[ip] += -gkb->zeta * gkb->d[ip];
    }

    /* Update error norm and test if one needs one more iteration */

    _gkb_cvg_test(gkb);

  }

  /* Return to the initial velocity formulation
   * u: = u_tilda + M^-1.(b_f + gamma.N^-1.b_c)
   * where M^-1.(b_f + gamma.N^-1.b_c) is stored in b_tilda */

# pragma omp parallel for if (gkb->n_u_dofs > CS_THR_MIN)
  for (cs_lnum_t iu = 0; iu < gkb->n_u_dofs; iu++)
    u_f[iu] = gkb->u_tilda[iu] + gkb->b_tilda[iu];

  int n_inner_iter = gkb->algo->n_inner_iter;

  /* Last step: Free temporary memory */

  _free_gkb_builder(&gkb);

  return  n_inner_iter;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Use the preconditioned Uzawa-CG algorithm to solve the saddle-point
 *         problem arising from CDO-Fb schemes for Stokes, Oseen and
 *         Navier-Stokes with a monolithic coupling
 *         This algorithm is based on Koko's paper "Uzawa conjugate gradient
 *         method for the Stokes problem: Matlab implementation with P1-iso-P2/
 *         P1 finite element"
 *
 * \param[in]      nsp      pointer to a cs_navsto_param_t structure
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in]      sh       pointer to a cs_cdo_system_helper_t structure
 * \param[in, out] slesp    pointer to a set of parameters to drive the SLES
 * \param[in, out] msles    pointer to a cs_cdofb_monolithic_sles_t structure
 *
 * \return the cumulated number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdofb_monolithic_uzawa_cg_solve(const cs_navsto_param_t       *nsp,
                                   const cs_equation_param_t     *eqp,
                                   const cs_cdo_system_helper_t  *sh,
                                   cs_param_sles_t               *slesp,
                                   cs_cdofb_monolithic_sles_t    *msles)
{
  int  _n_iter;

  assert(nsp != NULL && nsp->sles_param->strategy == CS_NAVSTO_SLES_UZAWA_CG);
  assert(sh != NULL && sh->n_blocks == 2);

  const cs_real_t  *B_op = msles->div_op;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_range_set_t  *range_set = cs_cdo_system_get_range_set(sh, 0);
  const cs_matrix_t  *matrix = cs_cdo_system_get_matrix(sh, 0);
  const cs_time_step_t  *ts = cs_glob_time_step;

  cs_real_t  *u_f = msles->u_f;
  cs_real_t  *p_c = msles->p_c;
  cs_real_t  *b_f = sh->rhs;
  cs_real_t  *b_c = sh->rhs + 3*msles->n_faces;

  /* Allocate and initialize the Uzawa-CG builder structure */

  cs_uza_builder_t  *uza = _init_uzawa_builder(nsp,
                                               0, /* grad-div scaling */
                                               3*msles->n_faces,
                                               msles->n_cells,
                                               cs_shared_quant);

  const cs_navsto_param_sles_t  *nslesp = nsp->sles_param;

  /* The Schur complement approximation (B.A^-1.Bt) is build and stored in the
     native format */

  cs_matrix_t  *smat = NULL;
  cs_real_t  *diag_smat = NULL, *xtra_smat = NULL;

  switch (nsp->sles_param->schur_approximation) {

  case CS_PARAM_SCHUR_DIAG_INVERSE:
    smat = _diag_schur_approximation(nsp, matrix, range_set,
                                     uza,
                                     &diag_smat, &xtra_smat);
    break;
  case CS_PARAM_SCHUR_LUMPED_INVERSE:
    smat = _invlumped_schur_approximation(nsp, eqp, matrix, range_set,
                                          slesp, msles, uza,
                                          &diag_smat, &xtra_smat);
    break;

  case CS_PARAM_SCHUR_MASS_SCALED:
    break; /* Nothing else to do. One relies on the pressure mass matrix */

  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid Schur approximation.",
              __func__);
  }

  cs_param_sles_t  *schur_slesp = nslesp->schur_sles_param;

  if (msles->schur_sles == NULL) /* has been defined by name */
    msles->schur_sles = cs_sles_find_or_add(-1, schur_slesp->name);

  /* Compute the first RHS: A.u0 = rhs = b_f - B^t.p_0 to solve */

  _apply_div_op_transpose(B_op, p_c, uza->rhs);

  if (range_set->ifs != NULL) {

    cs_interface_set_sum(range_set->ifs,
                         /* n_elts, stride, interlaced */
                         uza->n_u_dofs, 1, false, CS_REAL_TYPE,
                         uza->rhs);

    cs_interface_set_sum(range_set->ifs,
                         /* n_elts, stride, interlaced */
                         uza->n_u_dofs, 1, false, CS_REAL_TYPE,
                         b_f);

  }

# pragma omp parallel for if (uza->n_u_dofs > CS_THR_MIN)
  for (cs_lnum_t iu = 0; iu < uza->n_u_dofs; iu++)
    uza->rhs[iu] = b_f[iu] - uza->rhs[iu];

  /* Initial residual for rhs */

  double normalization = _get_fbvect_norm(uza->rhs);

  /* Compute the first velocity guess
   * Modify the tolerance in order to be more accurate on this step */

  char  *init_system_name = slesp->name;
  double  init_eps = slesp->eps;
  int  init_max_iter = slesp->n_max_iter;

  char  *system_name = NULL;
  BFT_MALLOC(system_name, strlen(eqp->name) + strlen(":init_guess") + 1, char);
  sprintf(system_name, "%s:init_guess", eqp->name);

  slesp->name = system_name;
  slesp->eps = fmin(slesp->eps, nslesp->il_algo_param.rtol);
  slesp->n_max_iter = CS_MAX(100, init_max_iter);

  cs_param_sles_update_cvg_settings(true, slesp); /* use the field id */

  _n_iter = cs_cdo_solve_scalar_system(uza->n_u_dofs,
                                       slesp,
                                       matrix,
                                       range_set,
                                       normalization,
                                       false, /* rhs_redux --> already done */
                                       msles->sles,
                                       u_f,
                                       uza->rhs);
  uza->algo->n_inner_iter += _n_iter;
  uza->algo->last_inner_iter += _n_iter;

  /* Set back the initial parameters */

  slesp->name = init_system_name;
  slesp->eps = init_eps;
  slesp->n_max_iter = init_max_iter;

  cs_param_sles_update_cvg_settings(true, slesp);

  /* Partial memory free */

  BFT_FREE(system_name);

  /* Set pointers used in this algorithm */

  cs_real_t  *gk = uza->gk;
  cs_real_t  *dk = uza->res_p;
  cs_real_t  *rk = uza->d__v;    /* P space */
  cs_real_t  *wk = uza->b_tilda; /* U space */
  cs_real_t  *dwk = uza->dzk;    /* P space */
  cs_real_t  *zk = uza->rhs;     /* P or U space */

  /* Compute the first residual rk0 (in fact the velocity divergence) */

  _apply_div_op(B_op, u_f, rk);

# pragma omp parallel for if (uza->n_p_dofs > CS_THR_MIN)
  for (cs_lnum_t ip = 0; ip < uza->n_p_dofs; ip++)
    rk[ip] = b_c[ip] - rk[ip];

  double div_l2_norm = _get_cbscal_norm(rk);

  /* Compute g0 as
   *   Solve smat.zk = r0
   *   g0 = alpha zk + nu Mp^-1 r0 */

  memset(zk, 0, sizeof(cs_real_t)*uza->n_p_dofs);

  _n_iter = 0;
  switch (nsp->sles_param->schur_approximation) {

  case CS_PARAM_SCHUR_DIAG_INVERSE:
  case CS_PARAM_SCHUR_LUMPED_INVERSE:
    _n_iter = cs_cdo_solve_scalar_cell_system(uza->n_p_dofs,
                                              schur_slesp,
                                              smat,
                                              div_l2_norm,
                                              msles->schur_sles,
                                              zk,
                                              rk);
    break;

  case CS_PARAM_SCHUR_MASS_SCALED:
    for (cs_lnum_t i = 0; i < msles->n_cells; i++)
      zk[i] = rk[i]/uza->inv_mp[i];
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid Schur approximation.",
              __func__);

  }

  uza->algo->n_inner_iter += _n_iter;
  uza->algo->last_inner_iter += _n_iter;

# pragma omp parallel for if (uza->n_p_dofs > CS_THR_MIN)
  for (cs_lnum_t ip = 0; ip < uza->n_p_dofs; ip++)
    gk[ip] = uza->alpha*zk[ip] + uza->inv_mp[ip]*rk[ip];

  /* dk0 <-- gk0 */

  memcpy(dk, gk, uza->n_p_dofs*sizeof(cs_real_t));

  uza->algo->normalization = cs_gdot(uza->n_p_dofs, rk, gk);
  uza->algo->res0 = uza->algo->normalization;
  uza->algo->res = uza->algo->res0;

  /* Main loop knowing g0, r0, d0, u0, p0 */
  /* ------------------------------------ */

  while (_uza_cg_cvg_test(uza)) {

    /* Sensitivity step: Compute wk as the solution of A.wk = B^t.dk */

    /* Define the rhs for this system */

    _apply_div_op_transpose(B_op, dk, uza->rhs);
    if (range_set->ifs != NULL)
      cs_interface_set_sum(range_set->ifs,
                           uza->n_u_dofs,
                           1, false, CS_REAL_TYPE, /* stride, interlaced */
                           uza->rhs);

    normalization = _get_fbvect_norm(uza->rhs);

    /* Solve A.wk = B^t.dk (should be -B^t this implies a sign modification
       during the update step) */

    memset(wk, 0, sizeof(cs_real_t)*uza->n_u_dofs);

    uza->algo->n_inner_iter
      += (uza->algo->last_inner_iter =
          cs_cdo_solve_scalar_system(uza->n_u_dofs,
                                     slesp,
                                     matrix,
                                     range_set,
                                     normalization,
                                     false, /* rhs_redux -->already done */
                                     msles->sles,
                                     wk,
                                     uza->rhs));

    _apply_div_op(B_op, wk, dwk); /* -B -w --> dwk has the right sign */

    normalization = _get_cbscal_norm(dwk);

    /* Solve smat.zk = dwk */

    memset(zk, 0, sizeof(cs_real_t)*uza->n_p_dofs);

    _n_iter = 0;
    switch (nsp->sles_param->schur_approximation) {

    case CS_PARAM_SCHUR_DIAG_INVERSE:
    case CS_PARAM_SCHUR_LUMPED_INVERSE:
      _n_iter = cs_cdo_solve_scalar_cell_system(uza->n_p_dofs,
                                                schur_slesp,
                                                smat,
                                                normalization,
                                                msles->schur_sles,
                                                zk,
                                                dwk);
      break;

    case CS_PARAM_SCHUR_MASS_SCALED:
      for (cs_lnum_t i = 0; i < msles->n_cells; i++)
        zk[i] = dwk[i]/uza->inv_mp[i];
      break;

    default:
      bft_error(__FILE__, __LINE__, 0, "%s: Invalid Schur approximation.",
                __func__);

    }

    uza->algo->n_inner_iter += _n_iter;
    uza->algo->last_inner_iter += _n_iter;

#   pragma omp parallel for if (uza->n_p_dofs > CS_THR_MIN)
    for (cs_lnum_t ip = 0; ip < uza->n_p_dofs; ip++)
      zk[ip] = uza->alpha*zk[ip] + uza->inv_mp[ip]*dwk[ip];

    /* Updates
     *  - Compute the rho_factor = <rk,dk> / <dk, dwk>
     *  - u(k+1) = u(k) + rho_factor * wk  --> --wk
     *  - p(k+1) = p(k) - rho_factor * dk
     *  - gk     = gk   - rho_factor * zk
     */

    double rho_factor_denum = cs_gdot(uza->n_p_dofs, dk, dwk);
    assert(fabs(rho_factor_denum) > 0);
    double  rho_factor = cs_gdot(uza->n_p_dofs, rk, dk) / rho_factor_denum;

#   pragma omp parallel for if (uza->n_u_dofs > CS_THR_MIN)
    for (cs_lnum_t iu = 0; iu < uza->n_u_dofs; iu++)
      u_f[iu] = u_f[iu] + rho_factor*wk[iu]; /* --wk */

#   pragma omp parallel for if (uza->n_p_dofs > CS_THR_MIN)
    for (cs_lnum_t ip = 0; ip < uza->n_p_dofs; ip++) {
      p_c[ip] = p_c[ip] - rho_factor*dk[ip];
      gk[ip]  = gk[ip]  - rho_factor*zk[ip];
      rk[ip]  = rk[ip]  - rho_factor*dwk[ip];
    }

    /* Conjugate gradient direction: update dk */

    double  beta_num = cs_gdot(uza->n_p_dofs, rk, gk);
    double  beta_factor = beta_num/uza->algo->res;

    uza->algo->res = beta_num;

    /* dk <-- gk + beta_factor * dk */

#   pragma omp parallel for if (uza->n_p_dofs > CS_THR_MIN)
    for (cs_lnum_t ip = 0; ip < uza->n_p_dofs; ip++)
      dk[ip] = gk[ip] + beta_factor*dk[ip];

  } /* End of main loop */

  /* Cumulated sum of iterations to solve the Schur complement and the velocity
     block (save before freeing the uzawa structure). */

  int n_inner_iter = uza->algo->n_inner_iter;

  /* Last step: Free temporary memory */

  BFT_FREE(diag_smat);
  BFT_FREE(xtra_smat);
  _free_uza_builder(&uza);

  return  n_inner_iter;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Use the Uzawa algorithm with an Augmented Lagrangian (ALU) technique
 *         in an incremental way to solve the saddle-point problem arising from
 *         CDO-Fb schemes for Stokes, Oseen and Navier-Stokes with a monolithic
 *         coupling
 *
 * \param[in]      nsp      pointer to a cs_navsto_param_t structure
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in]      sh       pointer to a cs_cdo_system_helper_t structure
 * \param[in, out] slesp    pointer to a set of parameters to drive the SLES
 * \param[in, out] msles    pointer to a cs_cdofb_monolithic_sles_t structure
 *
 * \return the cumulated number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdofb_monolithic_uzawa_al_incr_solve(const cs_navsto_param_t       *nsp,
                                        const cs_equation_param_t     *eqp,
                                        const cs_cdo_system_helper_t  *sh,
                                        cs_param_sles_t               *slesp,
                                        cs_cdofb_monolithic_sles_t    *msles)
{
  assert(nsp != NULL && nsp->sles_param->strategy == CS_NAVSTO_SLES_UZAWA_AL);

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_real_t  gamma = msles->graddiv_coef;
  const cs_real_t  *div_op = msles->div_op;
  const cs_range_set_t  *range_set = cs_cdo_system_get_range_set(sh, 0);
  const cs_matrix_t  *matrix = cs_cdo_system_get_matrix(sh, 0);

  cs_real_t  *u_f = msles->u_f;
  cs_real_t  *p_c = msles->p_c;
  cs_real_t  *b_f = sh->rhs;
  cs_real_t  *b_c = sh->rhs + 3*msles->n_faces;

  /* Allocate and initialize the ALU builder structure */

  cs_uza_builder_t  *uza = _init_uzawa_builder(nsp,
                                               gamma,
                                               3*msles->n_faces,
                                               msles->n_cells,
                                               quant);

  /* Transformation of the initial right-hand side */

  cs_real_t  *btilda_c = uza->d__v;
# pragma omp parallel for if (uza->n_p_dofs > CS_THR_MIN)
  for (cs_lnum_t ip = 0; ip < uza->n_p_dofs; ip++)
    btilda_c[ip] = uza->inv_mp[ip]*b_c[ip];

  _apply_div_op_transpose(div_op, btilda_c, uza->b_tilda);

  if (range_set->ifs != NULL) {

    cs_interface_set_sum(range_set->ifs,
                         /* n_elts, stride, interlaced */
                         uza->n_u_dofs, 1, false, CS_REAL_TYPE,
                         uza->b_tilda);

    cs_interface_set_sum(range_set->ifs,
                         /* n_elts, stride, interlaced */
                         uza->n_u_dofs, 1, false, CS_REAL_TYPE,
                         b_f);

  }

  /* Update the modify right-hand side: b_tilda = b_f + gamma*Dt.W^-1.b_c */

# pragma omp parallel for if (uza->n_u_dofs > CS_THR_MIN)
  for (cs_lnum_t iu = 0; iu < uza->n_u_dofs; iu++) {
    uza->b_tilda[iu] *= gamma;
    uza->b_tilda[iu] += b_f[iu];
  }

  /* Initialization */
  /* ============== */

  /* Compute the RHS for the Uzawa system: rhs = b_tilda - Dt.p_c */

  _apply_div_op_transpose(div_op, p_c, uza->rhs);

  if (range_set->ifs != NULL)
    cs_interface_set_sum(range_set->ifs,
                         /* n_elts, stride, interlaced */
                         uza->n_u_dofs, 1, false, CS_REAL_TYPE,
                         uza->rhs);

# pragma omp parallel for if (uza->n_u_dofs > CS_THR_MIN)
  for (cs_lnum_t iu = 0; iu < uza->n_u_dofs; iu++) {
    uza->rhs[iu] *= -1;
    uza->rhs[iu] += uza->b_tilda[iu];
  }

  /* Solve AL.u_f = rhs
   * Modifiy the tolerance in order to be more accurate on this step since
   * the accuracy at this step has an influence on the global accuracy
   */

  char  *init_system_name = slesp->name;
  double  init_eps = slesp->eps;

  char  *system_name = NULL;
  BFT_MALLOC(system_name, strlen(eqp->name) + strlen(":alu0") + 1, char);
  sprintf(system_name, "%s:alu0", eqp->name);

  slesp->name = system_name;
  slesp->eps = nsp->sles_param->il_algo_param.rtol;

  cs_real_t  normalization = cs_cdo_blas_square_norm_pfvp(uza->rhs);

  if (fabs(normalization) > 0)
    normalization = sqrt(normalization);
  else
    normalization = 1.0;

  uza->algo->n_inner_iter
    += (uza->algo->last_inner_iter =
        cs_cdo_solve_scalar_system(uza->n_u_dofs,
                                   slesp,
                                   matrix,
                                   range_set,
                                   normalization,
                                   false, /* rhs_redux */
                                   msles->sles,
                                   u_f,
                                   uza->rhs));

  /* Set back the initial parameters */

  slesp->name = init_system_name;
  slesp->eps = init_eps;

  /* Partial free */

  BFT_FREE(system_name);

  /* Main loop */
  /* ========= */

  cs_real_t  *delta_u = uza->b_tilda;
  cs_real_t  delta_u_l2 = normalization;

  /* Compute the divergence of u since this is a stopping criteria */

  _apply_div_op(div_op, u_f, uza->d__v);

  /* Update p_c = p_c - gamma * (D.u_f - b_c). Recall that B = -div
   * Compute the RHS for the Uzawa system: rhs = -gamma*B^T.W^-1.(B.u - g) */

# pragma omp parallel for if (uza->n_p_dofs > CS_THR_MIN)
  for (cs_lnum_t ip = 0; ip < uza->n_p_dofs; ip++) {
    uza->d__v[ip] -= b_c[ip];
    uza->res_p[ip] = uza->inv_mp[ip] * uza->d__v[ip];
    p_c[ip] += gamma * uza->res_p[ip];
  }

  while (_uza_incr_cvg_test(delta_u_l2, uza)) {

    /* Continue building the RHS */

    _apply_div_op_transpose(div_op, uza->res_p, uza->rhs);

    if (range_set->ifs != NULL)
      cs_interface_set_sum(range_set->ifs,
                           /* n_elts, stride, interlaced */
                           uza->n_u_dofs, 1, false, CS_REAL_TYPE,
                           uza->rhs);

#   pragma omp parallel for if (uza->n_u_dofs > CS_THR_MIN)
    for (cs_lnum_t iu = 0; iu < uza->n_u_dofs; iu++)
      uza->rhs[iu] *= -gamma;

    /* Solve AL.u_f = rhs */

    memset(delta_u, 0, sizeof(cs_real_t)*uza->n_u_dofs);

    uza->algo->n_inner_iter
      += (uza->algo->last_inner_iter =
          cs_cdo_solve_scalar_system(uza->n_u_dofs,
                                     slesp,
                                     matrix,
                                     range_set,
                                     delta_u_l2, /* normalization */
                                     false, /* rhs_redux */
                                     msles->sles,
                                     delta_u,
                                     uza->rhs));

    delta_u_l2 = cs_cdo_blas_square_norm_pfvp(delta_u);
    if (fabs(delta_u_l2) > 0)
      delta_u_l2 = sqrt(delta_u_l2);
    else
      delta_u_l2 = 1.0;

    /* Update the velocity */

#   pragma omp parallel for if (uza->n_u_dofs > CS_THR_MIN)
    for (cs_lnum_t iu = 0; iu < uza->n_u_dofs; iu++)
      u_f[iu] += delta_u[iu];

    /* Update the divergence */

    _apply_div_op(div_op, u_f, uza->d__v);

    /* Update p_c = p_c - gamma * (D.u_f - b_c). Recall that B = -div
     * Prepare the computation of the RHS for the Uzawa system:
     * rhs = -gamma*B^T.W^-1.(B.u - g) */

#   pragma omp parallel for if (uza->n_p_dofs > CS_THR_MIN)
    for (cs_lnum_t ip = 0; ip < uza->n_p_dofs; ip++) {
      uza->d__v[ip] -= b_c[ip];
      uza->res_p[ip] = uza->inv_mp[ip] * uza->d__v[ip];
      p_c[ip] += gamma * uza->res_p[ip];
    }

  } /* End of Uzawa iterations */

  int n_inner_iter = uza->algo->n_inner_iter;

  /* Last step: Free temporary memory */

  _free_uza_builder(&uza);

  return  n_inner_iter;
}

/*----------------------------------------------------------------------------*/

#undef _petsc_cmd

END_C_DECLS
