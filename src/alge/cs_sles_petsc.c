/*============================================================================
 * Sparse Linear Equation Solvers using PETSc
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PETSc headers
 *----------------------------------------------------------------------------*/

#include <petscconf.h>
#include <petscksp.h>
#include <petscmat.h>
#include <petscvec.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_log.h"
#include "cs_fp_exception.h"
#include "cs_halo.h"
#include "cs_matrix.h"
#include "cs_matrix_default.h"
#include "cs_matrix_petsc.h"
#include "cs_matrix_petsc_priv.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_sles.h"
#include "cs_sles_petsc.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_sles_petsc.c

  \brief handling of PETSc-based linear solvers

  \page sles_petsc PETSc-based linear solvers.

  \typedef cs_sles_petsc_setup_hook_t

  \brief Function pointer for user settings of a PETSc KSP solver setup.

  This function is called the end of the setup stage for a KSP solver.

  Note that using the advanced KSPSetPostSolve and KSPSetPreSolve functions,
  this also allows setting further function pointers for pre and post-solve
  operations (see the PETSc documentation).

  Note: if the context pointer is non-NULL, it must point to valid data
  when the selection function is called so that value or structure should
  not be temporary (i.e. local);

  \param[in, out]  context  pointer to optional (untyped) value or structure
  \param[in, out]  ksp      pointer to PETSc KSP context
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/* Basic per linear system options and logging */
/*---------------------------------------------*/

typedef struct _cs_sles_petsc_setup_t {

  KSP            ksp;                    /* Linear solver context */

  Mat            a;                      /* Linear system matrix */

  bool           share_a;                /* true if provided matrix already
                                            of PETSc type, false otherwise */

  double         r_norm;                 /* residual normalization */
  void          *cctx;                   /* convergence context */

} cs_sles_petsc_setup_t;

struct _cs_sles_petsc_t {

  /* Performance data */

  int                  n_setups;           /* Number of times system setup */
  int                  n_solves;           /* Number of times system solved */

  int                  n_iterations_last;  /* Number of iterations for last
                                              system resolution */
  int                  n_iterations_min;   /* Minimum number of iterations
                                              in system resolution history */
  int                  n_iterations_max;   /* Maximum number of iterations
                                              in system resolution history */
  int long long        n_iterations_tot;   /* Total accumulated number of
                                              iterations */

  cs_timer_counter_t   t_setup;            /* Total setup */
  cs_timer_counter_t   t_solve;            /* Total time used */

  /* Additional setup options */

  void                        *hook_context;   /* Optional user context */
  cs_sles_petsc_setup_hook_t  *setup_hook;     /* Post setup function */

  char     *matype_r;                      /* requested PETSc matrix type */
  MatType   matype;                        /* actual PETSc matrix type */

  /* Setup data */

  char                 *ksp_type;
  char                 *pc_type;
  KSPNormType           norm_type;

  cs_sles_petsc_setup_t   *setup_data;

};

/* Shell matrix context */
/*----------------------*/

typedef struct {

  const cs_matrix_t     *a;              /* Pointer to matrix */

  cs_matrix_row_info_t   r;              /* Access buffer */

} _mat_shell_t;

/*============================================================================
 *  Global variables
 *============================================================================*/

static int  _n_petsc_systems = 0;
static PetscViewer _viewer = NULL;
static PetscLogStage _log_stage[2];

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Export the linear system using PETSc Viewer mechanism
 *
 * parameters:
 *   name <-- name of the related linear system
 *   ksp  <-- Pointer to PETSc KSP structure
 *   rhs  <-- PETSc vector
 *----------------------------------------------------------------------------*/

static void
_export_petsc_system(const char   *name,
                     KSP           ksp,
                     Vec           b)
{
  const char *p = getenv("CS_PETSC_SYSTEM_VIEWER");

  if (p == NULL)
    return;

  /* Get system and preconditioner matrices */

  Mat a, pa;
  KSPGetOperators(ksp, &a, &pa);

  char  *filename = NULL;
  int len = strlen(name) + strlen("_matrix.dat") + 1;
  BFT_MALLOC(filename, len, char);

  if (strcmp(p, "BINARY") == 0) {

    PetscViewer viewer;

    /* Matrix */
    sprintf(filename, "%s_matrix.dat", name);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename, FILE_MODE_WRITE, &viewer);
    MatView(a, viewer);
    PetscViewerDestroy(&viewer);

    /* Right-hand side */
    sprintf(filename, "%s_rhs.dat", name);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename, FILE_MODE_WRITE, &viewer);
    VecView(b, viewer);
    PetscViewerDestroy(&viewer);

  }
  else if (strcmp(p, "ASCII") == 0) {

    PetscViewer viewer;

    /* Matrix */
    sprintf(filename, "%s_matrix.txt", name);
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer);
    MatView(a, viewer);
    PetscViewerDestroy(&viewer);

    /* Right-hand side */
    sprintf(filename, "%s_rhs.txt", name);
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer);
    VecView(b, viewer);
    PetscViewerDestroy(&viewer);

  }
  else if (strcmp(p, "MATLAB") == 0) {

    PetscViewer viewer;

    /* Matrix */
    sprintf(filename, "%s_matrix.m", name);
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer);
    PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
    MatView(a, viewer);
    PetscViewerPopFormat(viewer);
    PetscViewerDestroy(&viewer);

    /* Right-hand side */
    sprintf(filename, "%s_rhs.m", name);
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer);
    PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
    VecView(b, viewer);
    PetscViewerPopFormat(viewer);
    PetscViewerDestroy(&viewer);

  }

  BFT_FREE(filename);
}

/*----------------------------------------------------------------------------
 * Local matrix.vector product y = A.x with shell matrix.
 *
 * parameters:
 *   a  <-- Pointer to PETSc matrix structure
 *   x  <-- Multiplying vector values
 *   y  --> Resulting vector
 *----------------------------------------------------------------------------*/

static void
_shell_mat_mult(Mat  a,
                Vec  x,
                Vec  y)
{
  _mat_shell_t *sh;
  const PetscScalar *ax;
  PetscScalar *ay;

  assert(sizeof(PetscScalar) == sizeof(cs_real_t));

  MatShellGetContext(a, &sh);

  VecGetArrayRead(x, &ax);
  VecGetArray(y, &ay);

  cs_matrix_vector_multiply(sh->a, ax, ay);

  VecRestoreArray(x, &ax);
  VecRestoreArray(y, &ay);
}

/*----------------------------------------------------------------------------
 * Get shell matrix diagonal values
 *
 * parameters:
 *   a  <-- Pointer to PETSc matrix structure
 *   y  --> Resulting vector
 *----------------------------------------------------------------------------*/

static void
_shell_get_diag(Mat  a,
                Vec  y)
{
  _mat_shell_t *sh;
  PetscScalar *ay;

  assert(sizeof(PetscScalar) == sizeof(cs_real_t));

  MatShellGetContext(a, &sh);

  VecGetArray(y, &ay);

  cs_matrix_copy_diagonal(sh->a, ay);

  VecRestoreArray(y, &ay);
}

/*----------------------------------------------------------------------------
 * Get matrix row.
 *
 * This function assumes a CSR structure, but could be adapted quite
 * easily to MSR.
 *
 * parameters:
 *   a            <-- Pointer to PETSc matrix structure
 *   row          <-- Row to get
 *   nnz          --> Number of non-zeroes in the row
 *   cols         --> Column numbers
 *   vals         --> values
 *----------------------------------------------------------------------------*/

static void
_shell_get_row(Mat                 a,
               PetscInt            row,
               PetscInt           *nnz,
               const PetscInt     *cols[],
               const PetscScalar  *vals[])
{
  CS_UNUSED(nnz);

  _mat_shell_t *sh;

  MatShellGetContext(a, &sh);

  assert(sizeof(PetscScalar) == sizeof(cs_real_t));
  assert(sizeof(PetscInt) == sizeof(cs_lnum_t));

  cs_matrix_get_row(sh->a, row, &(sh->r));

  *cols = sh->r.col_id;
  *vals = sh->r.vals;
}

/*----------------------------------------------------------------------------
 * Duplicate matrix
 *
 * parameters:
 *   a  <-- Pointer to PETSc matrix structure
 *   op <-- Matrix duplication option
 *   m  --> Duplicate matrix
 *----------------------------------------------------------------------------*/

static void
_shell_mat_duplicate(Mat                  a,
                     MatDuplicateOption   op,
                     Mat                 *m)
{
  CS_UNUSED(op);

  _mat_shell_t *sh;

  MatShellGetContext(a, &sh);

  const cs_lnum_t n_rows = cs_matrix_get_n_rows(sh->a);

  /* Shell matrix */

  _mat_shell_t *shc;

  BFT_MALLOC(shc, 1, _mat_shell_t);
  shc->a = sh->a;
  cs_matrix_row_init(&(shc->r));

  MatCreateShell(PETSC_COMM_WORLD,
                 n_rows,
                 n_rows,
                 PETSC_DETERMINE,
                 PETSC_DETERMINE,
                 shc,
                 m);
}

/*----------------------------------------------------------------------------
 * Destroy matrix
 *
 * parameters:
 *   a  <-- Pointer to PETSc matrix structure
 *   op <-- Matrix duplication option
 *   m  --> Duplicate matrix
 *----------------------------------------------------------------------------*/

static void
_shell_mat_destroy(Mat                  a,
                   MatDuplicateOption   op,
                   Mat                 *m)
{
  CS_UNUSED(m);
  CS_UNUSED(op);

  _mat_shell_t *sh;

  MatShellGetContext(a, &sh);

  cs_matrix_row_finalize(&(sh->r));

  BFT_FREE(sh);
}

/*----------------------------------------------------------------------------
 * Convergence test using residual normalization.
 *
 * This test overloads KSPConvergedDefault, by changing the residual
 * normalization at the first time step.
 *
 * parameters:
 *   ksp    <-> KSP context
 *   n      <-- iteration id
 *   rnorm  <-- residual norm
 *   reason <-- convergence status
 *   ctx    <-> associated context
 *----------------------------------------------------------------------------*/

static PetscErrorCode
_cs_ksp_converged(KSP                  ksp,
                  PetscInt             n,
                  PetscReal            rnorm,
                  KSPConvergedReason  *reason,
                  void                *ctx)
{
  cs_sles_petsc_setup_t  *sd = ctx;

  if (!n)
    rnorm = sd->r_norm;

  return KSPConvergedDefault(ksp, n, rnorm, reason, sd->cctx);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * User function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Function pointer for user settings of a PETSc KSP solver setup.
 *
 * This function is called at the end of the setup stage for a KSP solver.
 *
 * Note that using the advanced KSPSetPostSolve and KSPSetPreSolve functions,
 * this also allows setting further function pointers for pre and post-solve
 * operations (see the PETSc documentation).
 *
 * Note: if the context pointer is non-NULL, it must point to valid data
 * when the selection function is called so that value or structure should
 * not be temporary (i.e. local);
 *
 * parameters:
 *   context <-> pointer to optional (untyped) value or structure
 *   ksp     <-> pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

void
cs_user_sles_petsc_hook(void               *context,
                        void               *ksp)
{
  CS_UNUSED(context);
  CS_UNUSED(ksp);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize PETSc if needed (calls cs_matrix_petsc_ensure_init).
 *----------------------------------------------------------------------------*/

void
cs_sles_petsc_init(void)
{
  cs_matrix_petsc_ensure_init();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define and associate a PETSc linear system solver
 *        for a given field or equation name.
 *
 * If this system did not previously exist, it is added to the list of
 * "known" systems. Otherwise, its definition is replaced by the one
 * defined here.
 *
 * This is a utility function: if finer control is needed, see
 * \ref cs_sles_define and \ref cs_sles_petsc_create.
 *
 * In case of rotational periodicity for a block (non-scalar) matrix,
 * the matrix type will be forced to MATSHELL ("shell") regardless
 * of the option used.
 *
 * Note that this function returns a pointer directly to the iterative solver
 * management structure. This may be used to set further options.
 * If needed, \ref cs_sles_find may be used to obtain a pointer to the matching
 * \ref cs_sles_t container.
 *
 * \param[in]      f_id          associated field id, or < 0
 * \param[in]      name          associated name if f_id < 0, or NULL
 * \param[in]      matrix_type   PETSc matrix type
 * \param[in]      setup_hook    pointer to optional setup epilogue function
 * \param[in,out]  context       pointer to optional (untyped) value or
 *                               structure for setup_hook, or NULL
 *
 * \return  pointer to newly created iterative solver info object.
 */
/*----------------------------------------------------------------------------*/

cs_sles_petsc_t *
cs_sles_petsc_define(int                          f_id,
                     const char                  *name,
                     const char                  *matrix_type,
                     cs_sles_petsc_setup_hook_t  *setup_hook,
                     void                        *context)
{
  cs_sles_petsc_t * c = cs_sles_petsc_create(matrix_type,
                                             setup_hook,
                                             context);

  cs_sles_t *sc = cs_sles_define(f_id,
                                 name,
                                 c,
                                 "cs_sles_petsc_t",
                                 cs_sles_petsc_setup,
                                 cs_sles_petsc_solve,
                                 cs_sles_petsc_free,
                                 cs_sles_petsc_log,
                                 cs_sles_petsc_copy,
                                 cs_sles_petsc_destroy);

  cs_sles_set_error_handler(sc,
                            cs_sles_petsc_error_post_and_abort);

  return c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create PETSc linear system solver info and context.
 *
 * In case of rotational periodicity for a block (non-scalar) matrix,
 * the matrix type will be forced to MATSHELL ("shell") regardless
 * of the option used.
 *
 * \param[in]      matrix_type   PETSc matrix type
 * \param[in]      setup_hook    pointer to optional setup epilogue function
 * \param[in,out]  context       pointer to optional (untyped) value or
 *                               structure for setup_hook, or NULL
 *
 * \return  pointer to associated linear system object.
 */
/*----------------------------------------------------------------------------*/

cs_sles_petsc_t *
cs_sles_petsc_create(const char                  *matrix_type,
                     cs_sles_petsc_setup_hook_t  *setup_hook,
                     void                        *context)

{
  cs_sles_petsc_t *c;

  PetscBool is_initialized;

  PetscInitialized(&is_initialized);
  if (is_initialized == PETSC_FALSE) {
#if defined(HAVE_MPI)
    if (cs_glob_n_ranks > 1)
      PETSC_COMM_WORLD = cs_glob_mpi_comm;
    else
      PETSC_COMM_WORLD = MPI_COMM_SELF;
#endif
    PetscInitializeNoArguments();
  }

  if (_viewer == NULL) {
    PetscLogStageRegister("Linear system setup", _log_stage);
    PetscLogStageRegister("Linear system solve", _log_stage + 1);
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, "petsc.log", &_viewer);
#if PETSC_VERSION_GE(3,7,0)
    PetscLogDefaultBegin();
#else
    PetscLogBegin();
#endif
  }

  _n_petsc_systems += 1;

  BFT_MALLOC(c, 1, cs_sles_petsc_t);
  c->n_setups = 0;
  c->n_solves = 0;
  c->n_iterations_last = 0;
  c->n_iterations_min = 0;
  c->n_iterations_max = 0;
  c->n_iterations_tot = 0;

  CS_TIMER_COUNTER_INIT(c->t_setup);
  CS_TIMER_COUNTER_INIT(c->t_solve);

  /* Options */

  c->hook_context = context;
  c->setup_hook = setup_hook;

  /* Setup data */

  PetscStrallocpy(matrix_type, &(c->matype_r));
  c->matype = NULL;

  c->setup_data = NULL;

  c->ksp_type = NULL;
  c->pc_type = NULL;
  c->norm_type = KSP_NORM_DEFAULT;

  return c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create PETSc linear system solver info and context
 *        based on existing info and context.
 *
 * \param[in]  context  pointer to reference info and context
 *                     (actual type: cs_sles_petsc_t  *)
 *
 * \return  pointer to newly created solver info object.
 *          (actual type: cs_sles_petsc_t  *)
 */
/*----------------------------------------------------------------------------*/

void *
cs_sles_petsc_copy(const void  *context)
{
  cs_sles_petsc_t *d = NULL;

  if (context != NULL) {
    const cs_sles_petsc_t *c = context;
    d = cs_sles_petsc_create(c->matype_r,
                             c->setup_hook,
                             c->hook_context);
  }

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy PETSc linear system solver info and context.
 *
 * \param[in, out]  context  pointer to iterative solver info and context
 *                           (actual type: cs_sles_petsc_t  **)
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_petsc_destroy(void  **context)
{
  cs_sles_petsc_t *c = (cs_sles_petsc_t *)(*context);
  if (c != NULL) {

    /* Free local strings */

    if (c->matype_r != NULL)
      PetscFree(c->matype_r);
    if (c->matype != NULL)
      PetscFree(c->matype);

    if (c->ksp_type != NULL)
      PetscFree(c->ksp_type);
    if (c->pc_type != NULL)
      PetscFree(c->pc_type);

    /* Free structure */

    cs_sles_petsc_free(c);
    BFT_FREE(c);
    *context = c;

    _n_petsc_systems -= 1;
    if (_n_petsc_systems == 0) {
      PetscLogView(_viewer);
      PetscViewerDestroy(&_viewer);
      cs_matrix_petsc_finalize();
    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Setup PETSc linear equation solver.
 *
 * \param[in, out]  context    pointer to iterative solver info and context
 *                             (actual type: cs_sles_petsc_t  *)
 * \param[in]       name       pointer to system name
 * \param[in]       a          associated matrix
 * \param[in]       verbosity  associated verbosity
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_petsc_setup(void               *context,
                    const char         *name,
                    const cs_matrix_t  *a,
                    int                 verbosity)
{
  cs_timer_t t0;
  t0 = cs_timer_time();

  PetscLogStagePush(_log_stage[0]);

  cs_sles_petsc_t  *c = context;
  cs_sles_petsc_setup_t *sd = c->setup_data;

  if (sd == NULL) {
    BFT_MALLOC(c->setup_data, 1, cs_sles_petsc_setup_t);
    sd = c->setup_data;
  }

  const cs_matrix_type_t cs_mat_type = cs_matrix_get_type(a);
  const PetscInt n_rows = cs_matrix_get_n_rows(a);
  const PetscInt db_size = cs_matrix_get_diag_block_size(a);
  const PetscInt eb_size = cs_matrix_get_extra_diag_block_size(a);
  const cs_halo_t *halo = cs_matrix_get_halo(a);

  bool have_perio = false;
  if (halo != NULL) {
    if (halo->n_transforms > 0)
      have_perio = true;
  }

  sd->share_a = false;

  /* Check if the matrix is already a PETSc matrix */

  if (strncmp(cs_matrix_get_type_name(a), "PETSc", 5) == 0) {
    cs_matrix_coeffs_petsc_t *coeffs = cs_matrix_petsc_get_coeffs(a);
    sd->a = coeffs->hm;
    sd->share_a = true;
  }

  /* Shell matrix */

  else if (   strcmp(c->matype_r, MATSHELL) == 0
           || (have_perio && n_rows > 1)
           || cs_mat_type == CS_MATRIX_NATIVE) {

    _mat_shell_t *sh;

    BFT_MALLOC(sh, 1, _mat_shell_t);
    sh->a = a;
    cs_matrix_row_init(&(sh->r));

    MatCreateShell(PETSC_COMM_WORLD,
                   n_rows*db_size,
                   n_rows*db_size,
                   PETSC_DECIDE,
                   PETSC_DECIDE,
                   sh,
                   &(sd->a));

    MatShellSetOperation(sd->a, MATOP_MULT,
                         (void(*)(void))_shell_mat_mult);
    MatShellSetOperation(sd->a, MATOP_GET_DIAGONAL,
                         (void(*)(void))_shell_get_diag);
    MatShellSetOperation(sd->a, MATOP_GET_ROW,
                         (void(*)(void))_shell_get_row);
    MatShellSetOperation(sd->a, MATOP_DUPLICATE,
                         (void(*)(void))_shell_mat_duplicate);
    MatShellSetOperation(sd->a, MATOP_DESTROY,
                         (void(*)(void))_shell_mat_destroy);

  }
  else if (   db_size == 1 && cs_mat_type == CS_MATRIX_CSR
           && (   strcmp(c->matype_r, MATMPIAIJ) == 0
               || (   strcmp(c->matype_r, MATAIJ) == 0
                   && cs_glob_n_ranks > 1))) {

    const cs_gnum_t *grow_id = cs_matrix_get_block_row_g_id(a);

    PetscInt *col_gid;
    const cs_lnum_t *a_row_index, *a_col_id;
    const cs_real_t *a_val;

    cs_matrix_get_csr_arrays(a, &a_row_index, &a_col_id, &a_val);

    BFT_MALLOC(col_gid, a_row_index[n_rows], PetscInt);

    for (cs_lnum_t j = 0; j < n_rows; j++) {
      for (cs_lnum_t i = a_row_index[j]; i < a_row_index[j+1]; ++i)
        col_gid[i] = grow_id[a_col_id[i]];
    }

    PetscInt     *row_index = a_row_index;
    PetscScalar  *val = a_val;

    if (sizeof(PetscInt) != sizeof(cs_lnum_t)) {
      BFT_MALLOC(row_index, n_rows, PetscInt);
      for (cs_lnum_t i = 0; i < n_rows; i++)
        row_index[i] = a_row_index[i];
    }

    if (sizeof(PetscScalar) != sizeof(cs_real_t)) {
      const cs_lnum_t val_size = a_row_index[n_rows];
      BFT_MALLOC(val, val_size, PetscScalar);
      for (cs_lnum_t i = 0; i < val_size; i++)
        val[i] = a_val[i];
    }

    /* Matrix */

    MatCreateMPIAIJWithArrays(PETSC_COMM_WORLD,  /* Communicator */
                              n_rows,            /* Number of local rows */
                              n_rows,            /* Number of local columns */
                              PETSC_DECIDE,      /* Number of global rows */
                              PETSC_DECIDE,      /* Number of global columns */
                              row_index,         /* Row indices */
                              col_gid,           /* Column indices */
                              val,               /* Matrix value */
                              &(sd->a));         /* Petsc Matrix */

    if (sizeof(PetscScalar) != sizeof(cs_real_t))
      BFT_FREE(val);

    if (sizeof(PetscInt) != sizeof(cs_lnum_t))
      BFT_FREE(row_index);

    BFT_FREE(col_gid);

  }
  else if (   sizeof(PetscInt) == sizeof(cs_lnum_t)
           && db_size == 1 && cs_mat_type == CS_MATRIX_CSR
           && (   strcmp(c->matype_r, MATSEQAIJ) == 0
               || (   strcmp(c->matype_r, MATAIJ) == 0
                   && cs_glob_n_ranks == 1))) {

    const cs_lnum_t *a_row_index, *a_col_id;
    const cs_real_t *a_val;

    cs_matrix_get_csr_arrays(a, &a_row_index, &a_col_id, &a_val);

    /* Matrix */

    MatCreateSeqAIJWithArrays(PETSC_COMM_SELF,  /* Communicator */
                              n_rows,           /* Number of local rows */
                              n_rows,           /* Number of local columns */
                              a_row_index,      /* Row indices */
                              a_col_id,         /* Column indices */
                              a_val,            /* Matrix value */
                              &(sd->a));        /* Petsc Matrix */

  }

  else {

    assert(cs_mat_type != CS_MATRIX_NATIVE);

    const cs_gnum_t *grow_id = cs_matrix_get_block_row_g_id(a);

    MatCreate(PETSC_COMM_WORLD, &(sd->a));
    MatSetType(sd->a, c->matype_r);
    MatSetSizes(sd->a,
                n_rows,            /* Number of local rows */
                n_rows,            /* Number of local columns */
                PETSC_DECIDE,      /* Number of global rows */
                PETSC_DECIDE);     /* Number of global columns */

    /* Preallocate */

    PetscInt *d_nnz, *o_nnz;
    BFT_MALLOC(d_nnz, n_rows*db_size, PetscInt);
    BFT_MALLOC(o_nnz, n_rows*db_size, PetscInt);

    if (cs_mat_type == CS_MATRIX_CSR || cs_mat_type == CS_MATRIX_MSR) {

      const cs_lnum_t *a_row_index, *a_col_id;
      const cs_real_t *a_val;
      const cs_real_t *d_val = NULL;

      if (cs_mat_type == CS_MATRIX_CSR) {

        cs_matrix_get_csr_arrays(a, &a_row_index, &a_col_id, &a_val);

        for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++) {
          for (cs_lnum_t kk = 0; kk < db_size; kk++) {
            d_nnz[row_id*db_size + kk] = 0;
            o_nnz[row_id*db_size + kk] = 0;
          }
        }

      }
      else {

        cs_matrix_get_msr_arrays(a, &a_row_index, &a_col_id, &d_val, &a_val);

        for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++) {
          for (cs_lnum_t kk = 0; kk < db_size; kk++) {
            d_nnz[row_id*db_size + kk] = db_size;
            o_nnz[row_id*db_size + kk] = 0;
          }
        }

      }

      for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++) {
        for (cs_lnum_t i = a_row_index[row_id]; i < a_row_index[row_id+1]; i++) {
          if (a_col_id[i] < n_rows) {
            for (cs_lnum_t kk = 0; kk < db_size; kk++)
              d_nnz[row_id*db_size + kk] += eb_size;
          }
          else {
            for (cs_lnum_t kk = 0; kk < db_size; kk++)
              o_nnz[row_id*db_size + kk] += eb_size;
          }
        }
      }

    }

    else
      bft_error(__FILE__, __LINE__, 0,
                _("Matrix type %s with block size %d for system \"%s\"\n"
                  "is not usable by PETSc."),
                cs_matrix_get_type_name(a), (int)db_size, name);

    /* Now preallocate matrix */

    MatSeqAIJSetPreallocation(sd->a, 0, d_nnz);
    MatMPIAIJSetPreallocation(sd->a, 0, d_nnz, 0, o_nnz);

    BFT_FREE(o_nnz);
    BFT_FREE(d_nnz);

    /* Now set matrix values, depending on type */

    if (cs_mat_type == CS_MATRIX_CSR || cs_mat_type == CS_MATRIX_MSR) {

      const cs_lnum_t *a_row_index, *a_col_id;
      const cs_real_t *a_val;
      const cs_real_t *d_val = NULL;

      PetscInt m = 1, n = 1;

      if (cs_mat_type == CS_MATRIX_CSR)
        cs_matrix_get_csr_arrays(a, &a_row_index, &a_col_id, &a_val);

      else {

        cs_matrix_get_msr_arrays(a, &a_row_index, &a_col_id, &d_val, &a_val);

        const cs_lnum_t b_size = cs_matrix_get_diag_block_size(a);
        const cs_lnum_t b_size_2 = b_size * b_size;

        for (cs_lnum_t b_id = 0; b_id < n_rows; b_id++) {
          for (cs_lnum_t ii = 0; ii < db_size; ii++) {
            for (cs_lnum_t jj = 0; jj < db_size; jj++) {
              PetscInt idxm[] = {grow_id[b_id*db_size + ii]};
              PetscInt idxn[] = {grow_id[b_id*db_size + jj]};
              PetscScalar v[] = {d_val[b_id*b_size_2 + ii*b_size + jj]};
              MatSetValues(sd->a, m, idxm, n, idxn, v, INSERT_VALUES);
            }
          }
        }

      }

      const cs_lnum_t b_size = cs_matrix_get_extra_diag_block_size(a);
      const cs_lnum_t b_size_2 = b_size * b_size;

      if (b_size == 1) {

        for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++) {
          for (cs_lnum_t i = a_row_index[row_id]; i < a_row_index[row_id+1]; i++) {
            cs_lnum_t c_id = a_col_id[i];
            for (cs_lnum_t kk = 0; kk < db_size; kk++) {
              PetscInt idxm[] = {grow_id[row_id*db_size + kk]};
              PetscInt idxn[] = {grow_id[c_id*db_size + kk]};
              PetscScalar v[] = {a_val[i]};
              MatSetValues(sd->a, m, idxm, n, idxn, v, INSERT_VALUES);
            }
          }
        }

      }
      else {

        for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++) {
          for (cs_lnum_t i = a_row_index[row_id]; i < a_row_index[row_id+1]; i++) {
            cs_lnum_t c_id = a_col_id[i];
            for (cs_lnum_t ii = 0; ii < db_size; ii++) {
              PetscInt idxm[] = {grow_id[row_id*db_size + ii]};
              for (cs_lnum_t jj = 0; jj < db_size; jj++) {
                PetscInt idxn[] = {grow_id[c_id*db_size + jj]};
                PetscScalar v[] = {d_val[i*b_size_2 + ii*b_size + jj]};
                MatSetValues(sd->a, m, idxm, n, idxn, v, INSERT_VALUES);
              }
            }
          }
        }

      }

    }

    MatAssemblyBegin(sd->a, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(sd->a, MAT_FINAL_ASSEMBLY);

  }

  /* Solver */

  KSPCreate(PETSC_COMM_WORLD, &(sd->ksp));
  KSPSetFromOptions(sd->ksp);
  KSPSetOperators(sd->ksp, sd->a, sd->a);

  KSPSetTolerances(sd->ksp, PETSC_DEFAULT, 1e-30, 1.e10, 10000);
  sd->r_norm = -1;

  KSPConvergedDefaultCreate(&(sd->cctx));
  KSPSetConvergenceTest(sd->ksp,
                        _cs_ksp_converged,
                        sd,
                        NULL);

  if (c->setup_hook != NULL)
    c->setup_hook(c->hook_context, sd->ksp);

  /* KSPSetup could be called here for better separation of setup/solve
     logging, but calling it systematically seems to cause issues
     at least with the performance of the GAMG preconditioner
     (possibly calling unneeded operations). So we avoid it for now,
     noting that the user always has to option of calling it at the
     end of the setup hook. */

  /* KSPSetUp(sd->ksp); */

  if (verbosity > 0)
    KSPView(sd->ksp, PETSC_VIEWER_STDOUT_WORLD);

  if (c->matype == NULL) {
    MatType matype;
    MatGetType(sd->a, &matype);
    PetscStrallocpy(matype, &(c->matype));
  }

  if (c->ksp_type == NULL) {
    KSPType ksptype;
    KSPGetType(sd->ksp, &ksptype);
    PetscStrallocpy(ksptype, &c->ksp_type);
    KSPGetNormType(sd->ksp, &(c->norm_type));
  }

  if (c->pc_type == NULL) {
    PC  pc;
    PCType pctype;
    KSPGetPC(sd->ksp, &pc);
    PCGetType(pc, &pctype);
    PetscStrallocpy(pctype, &c->pc_type);
  }

  PetscLogStagePop();

  /* Update return values */
  c->n_setups += 1;

  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(c->t_setup), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Call PETSc linear equation solver.
 *
 * \param[in, out]  context        pointer to iterative solver info and context
 *                                 (actual type: cs_sles_petsc_t  *)
 * \param[in]       name           pointer to system name
 * \param[in]       a              matrix
 * \param[in]       verbosity      associated verbosity
 * \param[in]       precision      solver precision
 * \param[in]       r_norm         residual normalization
 * \param[out]      n_iter         number of "equivalent" iterations
 * \param[out]      residual       residual
 * \param[in]       rhs            right hand side
 * \param[in, out]  vx             system solution
 * \param[in]       aux_size       number of elements in aux_vectors (in bytes)
 * \param           aux_vectors    optional working area
 *                                 (internal allocation if NULL)
 *
 * \return  convergence state
 */
/*----------------------------------------------------------------------------*/

cs_sles_convergence_state_t
cs_sles_petsc_solve(void                *context,
                    const char          *name,
                    const cs_matrix_t   *a,
                    int                  verbosity,
                    double               precision,
                    double               r_norm,
                    int                 *n_iter,
                    double              *residual,
                    const cs_real_t     *rhs,
                    cs_real_t           *vx,
                    size_t               aux_size,
                    void                *aux_vectors)
{
  CS_UNUSED(aux_size);
  CS_UNUSED(aux_vectors);

  cs_sles_convergence_state_t cvg = CS_SLES_ITERATING;

  cs_timer_t t0;
  t0 = cs_timer_time();

  cs_sles_petsc_t  *c = context;
  cs_sles_petsc_setup_t  *sd = c->setup_data;

  if (sd == NULL) {
    cs_sles_petsc_setup(c, name, a, verbosity);
    sd = c->setup_data;
  }

  PetscReal rtol, abstol, dtol;
  PetscInt  maxits;
  Vec       x, b;

  PetscLogStagePush(_log_stage[0]);

  KSPGetTolerances(sd->ksp, &rtol, &abstol, &dtol, &maxits);
  KSPSetTolerances(sd->ksp, precision, abstol, dtol, maxits);
  sd->r_norm = r_norm;

  PetscLogStagePop();
  PetscLogStagePush(_log_stage[1]);

  PetscInt       its;
  PetscScalar    _residual;
  const cs_lnum_t n_rows = cs_matrix_get_n_rows(a);
  const cs_lnum_t n_cols = cs_matrix_get_n_columns(a);
  const cs_lnum_t db_size = cs_matrix_get_diag_block_size(a);

  const PetscInt _n_rows = n_rows*db_size;
  cs_matrix_coeffs_petsc_t *coeffs = NULL;

  if (sd->share_a) {

    coeffs = cs_matrix_petsc_get_coeffs(a);

    x = coeffs->hx;
    b = coeffs->hy;

    PetscReal *_x = NULL, *_b = NULL;

    VecGetArray(x, &_x);
    VecGetArray(b, &_b);

#   pragma omp parallel for if(_n_rows > CS_THR_MIN)
    for (PetscInt i = 0; i < _n_rows; i++) {
      _x[i] = vx[i];
      _b[i] = rhs[i];
    }

    VecRestoreArray(x, &_x);
    VecRestoreArray(b, &_b);

  }
  else {

    if (cs_glob_n_ranks > 1) {

      PetscInt nghost = (n_cols - n_rows)*db_size;
      PetscInt *ghosts;

      BFT_MALLOC(ghosts, nghost, PetscInt);

      for (PetscInt i = 0; i < nghost; i++)
        ghosts[i] = n_rows*db_size + i;

      /* Vector */

      VecCreateGhostWithArray(PETSC_COMM_WORLD,
                              n_rows*db_size,
                              PETSC_DECIDE,
                              nghost,
                              ghosts,
                              vx,
                              &x);

      VecCreateGhostWithArray(PETSC_COMM_WORLD,
                              n_rows*db_size,
                              PETSC_DECIDE,
                              nghost,
                              ghosts,
                              rhs,
                              &b);

      BFT_FREE(ghosts);

    }
    else {

      VecCreateSeqWithArray(PETSC_COMM_SELF,
                            1,
                            n_rows*db_size,
                            vx,
                            &x);

      VecCreateSeqWithArray(PETSC_COMM_SELF,
                            1,
                            n_rows*db_size,
                            rhs,
                            &b);

    }

  }

  /* Export the linear system with PETSc functions */

  if (getenv("CS_PETSC_SYSTEM_VIEWER") != NULL)
    _export_petsc_system(name, sd->ksp, b);

  /* Resolution */

  cs_fp_exception_disable_trap();

  KSPSolve(sd->ksp, b, x);

  cs_fp_exception_restore_trap();

  if (sd->share_a) {

    PetscReal *_x = NULL;
    VecGetArrayRead(x, &_x);

#   pragma omp parallel for if(_n_rows > CS_THR_MIN)
    for (PetscInt i = 0; i < _n_rows; i++) {
      _x[i] = vx[i];
    }

    VecRestoreArrayRead(x, &_x);

  }
  else {

    VecDestroy(&x);
    VecDestroy(&b);

  }

  if (verbosity > 0)
    KSPView(sd->ksp, PETSC_VIEWER_STDOUT_WORLD);

  KSPGetResidualNorm(sd->ksp, &_residual);
  KSPGetIterationNumber(sd->ksp, &its);

  KSPConvergedReason reason;
  KSPGetConvergedReason(sd->ksp, &reason);

  if (reason > KSP_CONVERGED_ITERATING)
    cvg = CS_SLES_CONVERGED;
  else if (reason < KSP_CONVERGED_ITERATING) {
    if (reason == KSP_DIVERGED_ITS)
      cvg = CS_SLES_MAX_ITERATION;
    else if (   reason == KSP_DIVERGED_BREAKDOWN
             || reason == KSP_DIVERGED_BREAKDOWN_BICG)
      cvg = CS_SLES_BREAKDOWN;
    else
      cvg = CS_SLES_DIVERGED;
  }

  *residual = _residual;
  *n_iter = its;

  /* Update return values */

  PetscLogStagePop();

  if (c->n_solves == 0)
    c->n_iterations_min = its;

  c->n_iterations_last = its;
  c->n_iterations_tot += its;
  if (c->n_iterations_min > its)
    c->n_iterations_min = its;
  if (c->n_iterations_max < its)
    c->n_iterations_max = its;
  c->n_solves += 1;
  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(c->t_solve), &t0, &t1);

  return cvg;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free PETSc linear equation solver setup context.
 *
 * This function frees resolution-related data, such as
 * buffers and preconditioning but does not free the whole context,
 * as info used for logging (especially performance data) is maintained.
 *
 * \param[in, out]  context  pointer to iterative solver info and context
 *                           (actual type: cs_sles_petsc_t  *)
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_petsc_free(void  *context)
{
  cs_timer_t t0;
  t0 = cs_timer_time();

  cs_sles_petsc_t  *c  = context;
  cs_sles_petsc_setup_t *sd = c->setup_data;

  if (sd != NULL) {

    PetscLogStagePush(_log_stage[0]);

    KSPGetNormType(sd->ksp, &(c->norm_type));
    KSPConvergedDefaultDestroy(sd->cctx);
    KSPDestroy(&(sd->ksp));

    if (sd->share_a == false)
      MatDestroy(&(sd->a));

    PetscLogStagePop();

  }

  if (c->setup_data != NULL)
    BFT_FREE(c->setup_data);

  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(c->t_setup), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Error handler for PETSc solver.
 *
 * In case of divergence or breakdown, this error handler outputs an error
 * message
 * It does nothing in case the maximum iteration count is reached.

 * \param[in, out]  sles           pointer to solver object
 * \param[in]       state          convergence state
 * \param[in]       a              matrix
 * \param[in]       rhs            right hand side
 * \param[in, out]  vx             system solution
 *
 * \return  false (do not attempt new solve)
 */
/*----------------------------------------------------------------------------*/

bool
cs_sles_petsc_error_post_and_abort(cs_sles_t                    *sles,
                                   cs_sles_convergence_state_t   state,
                                   const cs_matrix_t            *a,
                                   const cs_real_t              *rhs,
                                   cs_real_t                    *vx)
{
  CS_UNUSED(a);
  CS_UNUSED(rhs);
  CS_UNUSED(vx);

  if (state >= CS_SLES_BREAKDOWN)
    return false;

  const cs_sles_petsc_t  *c = cs_sles_get_context(sles);
  const char *name = cs_sles_get_name(sles);

  const char *error_type[] = {N_("divergence"), N_("breakdown")};
  int err_id = (state == CS_SLES_BREAKDOWN) ? 1 : 0;

  bft_error(__FILE__, __LINE__, 0,
            _("%s and %s preconditioner with PETSc: error (%s) solving for %s"),
            _(c->ksp_type),
            _(c->pc_type),
            _(error_type[err_id]),
            name);

  return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log sparse linear equation solver info.
 *
 * \param[in]  context   pointer to iterative solver info and context
 *                       (actual type: cs_sles_petsc_t  *)
 * \param[in]  log_type  log type
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_petsc_log(const void  *context,
                  cs_log_t     log_type)
{
  const cs_sles_petsc_t  *c = context;

  const char undef[] = N_("not instanciated");
  const char *s_type = (c->ksp_type != NULL) ? c->ksp_type : undef;
  const char *p_type = (c->pc_type != NULL) ? c->pc_type : undef;
  const char *m_type = (c->matype != NULL) ? c->matype : undef;
  char norm_type_name[32];

  switch(c->norm_type) {
  case KSP_NORM_NONE:
    strncpy(norm_type_name, "none", 31);
    break;
  case KSP_NORM_PRECONDITIONED:
    strncpy(norm_type_name, "preconditioned", 31);
    break;
  case KSP_NORM_UNPRECONDITIONED:
    strncpy(norm_type_name, "unpreconditioned", 31);
    break;
  case KSP_NORM_NATURAL:
    strncpy(norm_type_name, "natural", 31);
    break;
  default:
    snprintf(norm_type_name, 31, "%d", c->norm_type);
  }
  norm_type_name[31] = '\0';

  if (log_type == CS_LOG_SETUP) {

    cs_log_printf(log_type,
                  _("  Solver type:                       PETSc (%s)\n"
                    "    Preconditioning:                   %s\n"
                    "    Norm type:                         %s\n"
                    "    Matrix format:                     %s\n"),
                    s_type, p_type, norm_type_name, m_type);

  }
  else if (log_type == CS_LOG_PERFORMANCE) {

    int n_calls = c->n_solves;
    int n_it_min = c->n_iterations_min;
    int n_it_max = c->n_iterations_max;
    int n_it_mean = 0;

    if (n_calls > 0)
      n_it_mean = (int)(  c->n_iterations_tot
                        / ((unsigned long long)n_calls));

    cs_log_printf(log_type,
                  _("\n"
                    "  Solver type:                   PETSc (%s)\n"
                    "    Preconditioning:             %s\n"
                    "    Norm type:                   %s\n"
                    "    Matrix format:               %s\n"
                    "  Number of setups:              %12d\n"
                    "  Number of calls:               %12d\n"
                    "  Minimum number of iterations:  %12d\n"
                    "  Maximum number of iterations:  %12d\n"
                    "  Mean number of iterations:     %12d\n"
                    "  Total setup time:              %12.3f\n"
                    "  Total solution time:           %12.3f\n"),
                  s_type, p_type, norm_type_name, m_type,
                  c->n_setups, n_calls, n_it_min, n_it_max, n_it_mean,
                  c->t_setup.nsec*1e-9,
                  c->t_solve.nsec*1e-9);

  }
}

/*----------------------------------------------------------------------------
 * \brief Output the settings of a KSP structure
 *
 * \param[in]  ksp     Krylov SubSpace structure
 *----------------------------------------------------------------------------*/

void
cs_sles_petsc_log_setup(void  *ksp)
{
  PetscViewer  v;

  PetscViewerCreate(PETSC_COMM_WORLD, &v);
  PetscViewerSetType(v, PETSCVIEWERASCII);
  PetscViewerFileSetMode(v, FILE_MODE_APPEND);
  PetscViewerFileSetName(v, "petsc_setup.log");

  KSPView(ksp, v);
  PetscViewerDestroy(&v);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the parameters driving the termination of an iterative process
 *        associated to a KSP structure
 *
 * \param[in, out] context   pointer to iterative solver info and context
 *                           (actual type: cs_sles_petsc_t  *)
 * \param[in]      rtol      relative tolerance
 * \param[in]      atol      absolute tolerance
 * \param[in]      dtol      divergence tolerance
 * \param[in]      max_it    max. number of iterations
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_petsc_set_cvg_criteria(const void      *context,
                               double           rtol,
                               double           atol,
                               double           dtol,
                               int              max_it)
{
  const cs_sles_petsc_t  *c = context;
  if (c == NULL)
    return;

  cs_sles_petsc_setup_t   *sd = c->setup_data;
  if (sd == NULL)
    return; /* No need to continue. This will be done during the first call to
               the solve function */

  KSPSetTolerances(sd->ksp, rtol, atol, dtol, max_it);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return matrix type associated with PETSc linear system solver
 *        info and context.
 *
 * \param[in, out]  context  pointer to iterative solver info and context
 *                           (actual type: cs_sles_petsc_t  **)
 *
 * \return  pointer to matrix type name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_sles_petsc_get_mat_type(void  *context)
{
  const char *retval = NULL;
  cs_sles_petsc_t *c = (cs_sles_petsc_t *)context;
  if (c != NULL)
    retval = c->matype_r;

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print information on PETSc library.
 *
 * \param[in]  log_type  log type
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_petsc_library_info(cs_log_t  log_type)
{
  cs_log_printf(log_type,
                "    PETSc %d.%d.%d\n",
                PETSC_VERSION_MAJOR,
                PETSC_VERSION_MINOR,
                PETSC_VERSION_SUBMINOR);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
