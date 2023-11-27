/*============================================================================
 * Sparse Linear Equation Solver Preconditioner driver
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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"

#include "cs_base.h"
#include "cs_blas.h"
#include "cs_field.h"
#include "cs_log.h"
#include "cs_halo.h"
#include "cs_mesh.h"
#include "cs_matrix.h"
#include "cs_matrix_default.h"
#include "cs_matrix_util.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_sles_pc.h"
#include "cs_sles_pc_priv.h"

#if defined(HAVE_CUDA)
#include "cs_sles_pc_cuda.h"
#endif

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_sles_pc.c

  \brief Sparse linear equation solver preconditioner driver and
  simple preconditioners.

  As the options associated with different preconditioners may be very
  varied, this subsystem is based on the use of a series of callback
  functions which may be associated with a given preconditioner.

  To summarize, the functions here define a preconditioner
  driver, with the real work being done by functions bound to this model.
  The main intent is to help manage passing varied user options to the
  preconditioners, and handling consistency of logging.

  \enum cs_sles_pc_state_t

  \brief Convergence status indicator.

  \var CS_SLES_PC_DIVERGED
       The preconditioner has diverged
  \var CS_SLES_PC_BREAKDOWN
       The preconditioner has broken down, and cannot make any more progress
  \var CS_SLES_PC_MAX_ITERATION
       Maximum number of iterations has been reached, without reaching
       convergence
  \var CS_SLES_CONVERGED
       The preconditioner has converged

  \typedef  cs_sles_pc_get_type_t

  \brief  Function pointer returning the type name of a preconditioner
          context.

  \param[in]  context    pointer to solver context
  \param[in]  logging    if true, a name appropritate to logging
                         (possibly translated) is returned; if false,
                         a canonical name is returned.

  \typedef  cs_sles_pc_setup_t

  \brief  Function pointer for pre-resolution setup of a preconditioner
          context.

  This setup may include building a multigrid hierarchy, for example.

  \param[in, out]  context    pointer to solver context
  \param[in]       name       pointer to name of linear system
  \param[in]       a          matrix
  \param[in]       verbosity  associated verbosity

  \typedef  cs_sles_pc_tolerance_t

  \brief  Function pointer for setting of the required tolerance for
  preconditioners involving an iterative solver.

  This will usually not be relevant to non-iterative preconditioners,
  for which this type of function does not need to be defined.

  The preconditioner is considered to have converged when
  residual/r_norm <= precision, residual being the L2 norm of a.vx-rhs.

  \param[in, out]  context        pointer to solver context
  \param[in]       precision      solver precision
  \param[in]       r_norm         residual normalization

  \typedef  cs_sles_pc_apply_t

  \brief  Function pointer for application of a preconditioner.

  In cases where it is desired that the preconditioner modify a vector
  "in place", x_in should be set to NULL, and x_out contain the vector to
  be modified (\f$x_{out} \leftarrow M^{-1}x_{out})\f$).

  \param[in, out]  context        pointer to solver context
  \param[in]       x_in           input vector
  \param[in, out]  x_out          input/output vector

  \return  preconditioner convergence status

  \typedef  cs_sles_pc_free_t

  \brief  Function pointer for freeing of a preconditioner's context data.

  Note that this function should free resolution-related data, such as
  multigrid hierarchy and any other temporary arrays or
  objects required for application, but should not free the whole context,
  as info used for logging (especially performance data) should be
  maintained.

  \param[in, out]  context  pointer to solver context.

  \typedef  cs_sles_pc_log_t

  \brief  Function pointer for logging of preconditioner history
          and performance data.

  This function will indirectly  be called for each preconditioner when
  \ref cs_sles_finalize is called.

  \param[in]  context   pointer to preconditioner context
  \param[in]  log_type  log type

  \typedef  cs_sles_pc_clone_t

  \brief  Function pointer for creation of a preconditioner context based on
          the copy of another.

  The new context copies the settings of the copied context, but not
  its setup data and logged info, such as performance data.

  This type of function is optional, but enables associating different
  preconditioners to related systems (to differentiate logging) while using
  the same settings by default.

  \param[in]  context  context to clone

  \return  pointer to newly created context

  \typedef  cs_sles_pc_destroy_t

  Function pointer for destruction of a preconditioner context.

  This function should free all context data, and will be called for each
  system when \ref cs_sles_finalize is called.

  \param[in, out]  context  pointer to preconditioner context

  Simple preconditioners also provided include:

  - identity (null)
  - Jacobi
  - polynomial of degree 1
  - polynomial of degree 2

  Polynomial preconditioning is explained here:
  \a D being the diagonal part of matrix \a A and \a X its extra-diagonal
  part, it can be written \f$A=D(Id+D^{-1}X)\f$. Therefore
  \f$A^{-1}=(Id+D^{-1}X)^{-1}D^{-1}\f$. A series development of
  \f$Id+D^{-1}X\f$ can then be used which yields, symbolically,
  \f[
  Id+\sum\limits_{I=1}^{poly\_degree}\left(-D^{-1}X\right)^{I}
  \f]

  The Jacobi preconditioning is a polynomial preconditioning of degree 0.

  The efficiency of the polynomial preconditioning will vary depending
  on the system type. In most cases, Jacobi or degree 1 provide
  best results. Each polynomial preconditioning degree above 0 adds one
  matrix-vector product per inital matrix-vector product of the algorithm.
  Switching from diagonal to polynomial degree 1 often divides the number of
  required iterations by approximately 2, but each iteration then costs
  close to 2 times that of diagonal preconditoning (other vector operations
  are not doubled), so the net gain is often about 10%. Higher degree
  polynomials usually lead to diminishing returns, so to avoid the need
  for additional parameter setting functions, only degrees 1
  and 2 are provided here.

*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/* SIMD unit size to ensure SIMD alignement (2 to 4 required on most
 * current architectures, so 16 should be enough on most architectures) */

#define CS_SIMD_SIZE(s) (((s-1)/16+1)*16)

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/* Basic per linear system options and logging */
/*---------------------------------------------*/

struct _cs_sles_pc_t {

  void                     *context;         /* preconditioner context
                                                (options, state, logging) */

  cs_sles_pc_get_type_t    *get_type_func;   /* type string function */

  cs_sles_pc_setup_t       *setup_func;      /* solver setup function */
  cs_sles_pc_tolerance_t   *tolerance_func;  /* set tolerance function */
  cs_sles_pc_apply_t       *apply_func;      /* apply function */
  cs_sles_pc_free_t        *free_func;       /* free setup function */

  cs_sles_pc_log_t         *log_func;        /* logging function */

  cs_sles_pc_clone_t       *clone_func;      /* clone function */
  cs_sles_pc_destroy_t     *destroy_func;    /* destruction function */

};

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create a Polynomial preconditioner structure.
 *
 * returns:
 *   pointer to newly created preconditioner object.
 *----------------------------------------------------------------------------*/

static cs_sles_pc_poly_t *
_sles_pc_poly_create(void)
{
  cs_sles_pc_poly_t *pc;

  BFT_MALLOC(pc, 1, cs_sles_pc_poly_t);

#if defined(HAVE_ACCEL)
  pc->accelerated = false;
#endif

  pc->poly_degree = 0;

  pc->n_rows = 0;
  pc->n_cols = 0;
  pc->n_aux = 0;

  pc->ad_inv = NULL;
  pc->_ad_inv = NULL;

  pc->aux = NULL;

  return pc;
}

/*----------------------------------------------------------------------------
 * Function returning the type name of polynomial preconditioner context.
 *
 * parameters:
 *   context   <-- pointer to preconditioner context
 *   logging   <-- if true, logging description; if false, canonical name
 *----------------------------------------------------------------------------*/

static const char *
_sles_pc_poly_get_type(const void  *context,
                       bool         logging)
{
  const cs_sles_pc_poly_t  *c = context;

  assert(c->poly_degree > -2 && c->poly_degree < 3);

  if (logging == false) {
    static const char *t[] = {"none",
                              "jacobi",
                              "polynomial_degree_1",
                              "polynomial_degree_2"};
    return t[c->poly_degree + 1];
  }
  else {
    static const char *t[] = {N_("none"),
                              N_("Jacobi"),
                              N_("polynomial, degree 1"),
                              N_("polynomial, degree 2")};
    return _(t[c->poly_degree + 1]);
  }
}

/*----------------------------------------------------------------------------
 * Function for setup of a polynomial preconditioner context.
 *
 * parameters:
 *   context   <-> pointer to preconditioner context
 *   name      <-- pointer to name of associated linear system
 *   a         <-- matrix
 *   accel     <-- use accelerator version ?
 *   verbosity <-- associated verbosity
 *----------------------------------------------------------------------------*/

static void
_sles_pc_poly_setup(void               *context,
                    const char         *name,
                    const cs_matrix_t  *a,
                    bool                accel,
                    int                 verbosity)
{
  CS_UNUSED(name);
  CS_UNUSED(verbosity);

  cs_sles_pc_poly_t  *c = context;

#if defined(HAVE_ACCEL)
  c->accelerated = accel;
  cs_alloc_mode_t amode = (c->accelerated) ?
    CS_ALLOC_HOST_DEVICE_SHARED : CS_ALLOC_HOST;
#else
  CS_UNUSED(accel);
  cs_alloc_mode_t amode = CS_ALLOC_HOST;
#endif

  const cs_lnum_t db_size = cs_matrix_get_diag_block_size(a);

  cs_lnum_t n_rows_prev = c->n_rows;
  c->n_rows = cs_matrix_get_n_rows(a)*db_size;
  c->n_cols = cs_matrix_get_n_columns(a)*db_size;

  c->a = a;

  const cs_lnum_t n_rows = c->n_rows;

  if (c->n_rows > n_rows_prev) {
    CS_FREE_HD(c->_ad_inv);
    CS_MALLOC_HD(c->_ad_inv, n_rows, cs_real_t, amode);
  }
  c->ad_inv = c->_ad_inv;

  cs_matrix_copy_diagonal(a, c->_ad_inv);

# pragma omp parallel for if(n_rows > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_rows; i++)
    c->_ad_inv[i] = 1.0 / c->_ad_inv[i];


#if defined(HAVE_ACCEL)
  if (c->accelerated)
    cs_sync_h2d_future(c->_ad_inv);
#endif
}

/*----------------------------------------------------------------------------
 * Function for setup when poly_degree < 0
 *
 * parameters:
 *   context   <-> pointer to preconditioner context
 *   name      <-- pointer to name of associated linear system
 *   a         <-- matrix
 *   accel     <-- use accelerator version ?
 *   verbosity <-- associated verbosity
 *----------------------------------------------------------------------------*/

static void
_sles_pc_poly_setup_none(void               *context,
                         const char         *name,
                         const cs_matrix_t  *a,
                         bool                accel,
                         int                 verbosity)
{
  CS_UNUSED(name);
  CS_UNUSED(verbosity);

  cs_sles_pc_poly_t  *c = context;

#if defined(HAVE_ACCEL)
  c->accelerated = accel;
#else
  CS_UNUSED(accel);
#endif

  const cs_lnum_t db_size = cs_matrix_get_diag_block_size(a);

  c->n_rows = cs_matrix_get_n_rows(a)*db_size;
  c->n_cols = cs_matrix_get_n_columns(a)*db_size;

  c->a = a;

  c->ad_inv = NULL;
  c->_ad_inv = NULL;

  c->n_aux = 0;
  c->aux = NULL;
}

/*----------------------------------------------------------------------------
 * Function for application of a null-preconditioner.
 *
 * In cases where it is desired that the preconditioner modify a vector
 * "in place", x_in should be set to NULL, and x_out contain the vector to
 * be modified (\f$x_{out} \leftarrow M^{-1}x_{out})\f$).
 *
 * parameters:
 *   context       <-> pointer to preconditioner context
 *   x_in          <-- input vector
 *   x_out         <-> input/output vector
 *
 * returns:
 *   preconditioner application status
 *----------------------------------------------------------------------------*/

static cs_sles_pc_state_t
_sles_pc_poly_apply_none(void                *context,
                         const cs_real_t     *x_in,
                         cs_real_t           *x_out)
{
  if (x_in != NULL) {

    cs_sles_pc_poly_t  *c = context;
    const cs_lnum_t n_rows = c->n_rows;

#if defined(HAVE_CUDA)
    if (c->accelerated)
      return cs_sles_pc_cuda_apply_none(context, x_in, x_out);
#endif

#   pragma omp parallel for if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++)
      x_out[ii] = x_in[ii];
  }

  return CS_SLES_PC_CONVERGED;
}

/*----------------------------------------------------------------------------
 * Function for application of a Jacobi preconditioner.
 *
 * In cases where it is desired that the preconditioner modify a vector
 * "in place", x_in should be set to NULL, and x_out contain the vector to
 * be modified (\f$x_{out} \leftarrow M^{-1}x_{out})\f$).
 *
 * parameters:
 *   context       <-> pointer to preconditioner context
 *   x_in          <-- input vector
 *   x_out         <-> input/output vector
 *
 * returns:
 *   preconditioner application status
 *----------------------------------------------------------------------------*/

static cs_sles_pc_state_t
_sles_pc_poly_apply_jacobi(void                *context,
                           const cs_real_t     *x_in,
                           cs_real_t           *x_out)
{
  cs_sles_pc_poly_t  *c = context;

#if defined(HAVE_CUDA)
  if (c->accelerated)
    return cs_sles_pc_cuda_apply_jacobi(context, x_in, x_out);
#endif

  const cs_lnum_t n_rows = c->n_rows;
  const cs_real_t *restrict ad_inv = c->ad_inv;

  if (x_in != NULL) {
#   pragma omp parallel for if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++)
      x_out[ii] = x_in[ii] * ad_inv[ii];
  }
  else {
#   pragma omp parallel for if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++)
      x_out[ii] *= ad_inv[ii];
  }

  return CS_SLES_PC_CONVERGED;
}

/*----------------------------------------------------------------------------
 * Function for application of a polynomial preconditioner.
 *
 * In cases where it is desired that the preconditioner modify a vector
 * "in place", x_in should be set to NULL, and x_out contain the vector to
 * be modified (\f$x_{out} \leftarrow M^{-1}x_{out})\f$).
 *
 * parameters:
 *   context       <-> pointer to preconditioner context
 *   x_in          <-- input vector
 *   x_out         <-> input/output vector
 *
 * returns:
 *   preconditioner application status
 *----------------------------------------------------------------------------*/

static cs_sles_pc_state_t
_sles_pc_poly_apply_poly(void                *context,
                         const cs_real_t     *x_in,
                         cs_real_t           *x_out)
{
  cs_sles_pc_poly_t  *c = context;

#if defined(HAVE_CUDA)
  if (c->accelerated)
    return cs_sles_pc_cuda_apply_poly(context, x_in, x_out);
#endif

  const cs_lnum_t n_rows = c->n_rows;
  const cs_lnum_t n_aux = (x_in == NULL) ?
    CS_SIMD_SIZE(c->n_cols) + c->n_cols : c->n_cols;

  if (c->n_aux < n_aux) {
    c->n_aux = n_aux;
    cs_alloc_mode_t amode = cs_check_device_ptr(c->ad_inv);
    CS_FREE_HD(c->aux);
    CS_MALLOC_HD(c->aux, c->n_aux, cs_real_t, amode);
  }

  cs_real_t *restrict w = c->aux;
  const cs_real_t *restrict r = x_in;
  const cs_real_t *restrict ad_inv = c->ad_inv;

  if (x_in == NULL) {

    cs_real_t *restrict _r = c->aux + CS_SIMD_SIZE(c->n_cols);

#   pragma omp parallel for if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++)
      _r[ii] = x_out[ii];

    r = _r;

  }

# pragma omp parallel for if(n_rows > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < n_rows; ii++)
    x_out[ii] = r[ii] * ad_inv[ii];

  for (int deg_id = 1; deg_id <= c->poly_degree; deg_id++) {

    /* Compute Wk = (A-diag).Gk */

    cs_matrix_vector_multiply_partial(c->a, CS_MATRIX_SPMV_E, x_out, w);

#   pragma omp parallel for if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++)
      x_out[ii] = (r[ii] - w[ii]) * ad_inv[ii];

  }

  return CS_SLES_PC_CONVERGED;
}

/*----------------------------------------------------------------------------
 * Function for freeing of a polynomial preconditioner's context data.
 *
 * parameters:
 *   context <-> pointer to preconditioner context
 *----------------------------------------------------------------------------*/

static void
_sles_pc_poly_free(void  *context)
{
  cs_sles_pc_poly_t  *c = context;

  c->n_rows = 0;
  c->n_cols = 0;
  c->n_aux = 0;

  c->a = NULL;

  c->ad_inv = NULL;
  CS_FREE_HD(c->_ad_inv);
  CS_FREE_HD(c->aux);
}

/*----------------------------------------------------------------------------
 * Function for creation of a polynomial preconditioner context based on the
 * copy of another.
 *
 * The new context copies the settings of the copied context, but not
 * its setup data and logged info, such as performance data.
 *
 * This type of function is optional, but enables associating different
 * preconditioners to related systems (to differentiate logging) while using
 * the same settings by default.
 *
 * parameters:
 *   context  <-- context to clone
 *
 * returns:
 *   pointer to newly created context
 *----------------------------------------------------------------------------*/

static void *
_sles_pc_poly_clone(const void  *context)
{
  const cs_sles_pc_poly_t *c = (const cs_sles_pc_poly_t *)context;

  cs_sles_pc_poly_t *pc = _sles_pc_poly_create();

  pc->poly_degree = c->poly_degree;

  return pc;
}

/*----------------------------------------------------------------------------
 * Function pointer for destruction of a preconditioner context.
 *
 * This function should free all context data.
 *
 * parameters:
 *   context <-> pointer to preconditioner context
 *----------------------------------------------------------------------------*/

static void
_sles_pc_poly_destroy (void  **context)
{
  if (context != NULL) {
    _sles_pc_poly_free(*context);
    BFT_FREE(*context);
  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define sparse linear equation preconditioner.
 *
 * The context pointer is used to point to a structure adapted to the function
 * pointers given here, and combined with those functions, allows using
 * both built-in, external, or user-defined preconditioners.
 *
 * \param[in, out]  context         pointer to preconditioner context structure
 *                                  (cs_sles_pc subsystem becomes owner)
 * \param[in]       get_type_func   pointer to function returning type name
 * \param[in]       setup_func      pointer to preconditioner setup function
 * \param[in]       tolerance_func  pointer to tolerance setting functio
 * \param[in]       apply_func      pointer to preconditioner application
 *                                  function (also calls setup_func if not done
 *                                  yet or free_func called since last apply)
 * \param[in]       free_func       pointer function freeing system setup
 * \param[in]       log_func        pointer to system info logging function
                                    (optional, but recommended)
 * \param[in]       clone_func      pointer to settings clone function
 * \param[in]       destroy_func    pointer to function destroying
 *                                  preconditioner context
 *
 * \return  pointer to associated preconditioner object
 */
/*----------------------------------------------------------------------------*/

cs_sles_pc_t *
cs_sles_pc_define(void                    *context,
                  cs_sles_pc_get_type_t   *get_type_func,
                  cs_sles_pc_setup_t      *setup_func,
                  cs_sles_pc_tolerance_t  *tolerance_func,
                  cs_sles_pc_apply_t      *apply_func,
                  cs_sles_pc_free_t       *free_func,
                  cs_sles_pc_log_t        *log_func,
                  cs_sles_pc_clone_t      *clone_func,
                  cs_sles_pc_destroy_t    *destroy_func)
{
  cs_sles_pc_t  *pc;

  BFT_MALLOC(pc, 1, cs_sles_pc_t);

  /* Now define options */
  pc->context = context;
  pc->get_type_func = get_type_func;
  pc->setup_func = setup_func;
  pc->tolerance_func = tolerance_func;
  pc->apply_func = apply_func;
  pc->free_func = free_func;
  pc->log_func = log_func;
  pc->clone_func = clone_func;
  pc->destroy_func = destroy_func;

  return pc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy a sparse linear equation preconditioner.
 *
 *
 * \param[in, out]  pc   pointer to preconditioner context structure
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_pc_destroy(cs_sles_pc_t **pc)
{
  if (pc != NULL) {
    cs_sles_pc_t *_pc = *pc;
    if (_pc != NULL) {
      _pc->destroy_func(&(_pc->context));
      BFT_FREE(*pc);
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a new preconditioner context based on the copy of another.
 *
 * The intended use of this function is to allow associating different
 * preconditioners to related systems, so as to allow simultaneous setups
 * and differentiate logging, while using the same settings by default.
 *
 * If no preconditioner (i.e. NULL) is passed, it will return NULL.
 *
 * \param[in]       src   pointer to source preconditioner object
 *
 * \return  pointer to new preconditioner object, or NULL
 */
/*----------------------------------------------------------------------------*/

cs_sles_pc_t *
cs_sles_pc_clone(const cs_sles_pc_t  *src)
{
  if (src == NULL)
    return NULL;

  cs_sles_pc_t  *dest;

  BFT_MALLOC(dest, 1, cs_sles_pc_t);

  /* Now define options */

  dest->context = src->clone_func(src->context);
  dest->get_type_func = src->get_type_func;
  dest->setup_func = src->setup_func;
  dest->tolerance_func = src->tolerance_func;
  dest->apply_func = src->apply_func;
  dest->free_func = src->free_func;
  dest->log_func = src->log_func;
  dest->clone_func = src->clone_func;
  dest->destroy_func = src->destroy_func;

  return dest;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return type name of preconditioner context.
 *
 * The returned string is intended to help determine which type is associated
 * with the void * pointer returned by \ref cs_sles_pc_get_context for a given
 * preconditioner definition, so as to be able to call additional specific
 * functions beyond the generic functions assigned to a cs_sles_pc_t object.
 *
 * \param[in]  pc  pointer to preconditioner object
 *
 * \return  pointer to linear system preconditioner specific type name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_sles_pc_get_type(cs_sles_pc_t  *pc)
{
  if (pc != NULL)
    return pc->get_type_func(pc->context, false);
  else {
    static const char t[] = "none";
    return t;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return type name of preconditioner context.
 *
 * The returned string is intended mainly for logging.
 *
 * \param[in]  pc  pointer to preconditioner object
 */
/*----------------------------------------------------------------------------*/

const char *
cs_sles_pc_get_type_name(cs_sles_pc_t  *pc)
{
  if (pc != NULL)
    return pc->get_type_func(pc->context, true);
  else {
    static const char t[] = "none";
    return _(t);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointer to preconditioner context structure pointer.
 *
 * The context structure depends on the type of preconditioner used, which may
 * in turn be determined by the string returned by cs_sles_pc_get_type().
 * If may be used by appropriate functions specific to that type.
 *
 * \param[in]  pc  pointer to preconditioner object
 *
 * \return  pointer to preconditioner-specific info and context
 */
/*----------------------------------------------------------------------------*/

void *
cs_sles_pc_get_context(cs_sles_pc_t  *pc)
{
  void *c = NULL;

  if (pc != NULL)
    c = pc->context;

  return c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a pointer to the function used to apply a preconditioner.
 *
 * This allows calling the preconditioner with one less level of indirection.
 *
 * \param[in]  pc  pointer to preconditioner object
 *
 * \return  preconditioner apply function
 */
/*----------------------------------------------------------------------------*/

cs_sles_pc_apply_t *
cs_sles_pc_get_apply_func(const cs_sles_pc_t *pc)
{
  return pc->apply_func;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the required tolerance for preconditioners involving an
 *        iterative solver.
 *
 * This will usually not be relevant to non-iterative preconditioners,
 * in which case this is a no-op.
 *
 * If no options were previously provided for the matching system,
 * default options will be used.
 *
 * The system is considered to have converged when
 * residual/r_norm <= precision, residual being the L2 norm of a.vx-rhs.
 *
 * \param[in, out]  pc             pointer to preconditioner object
 * \param[in]       precision      preconditioner precision
 * \param[in]       r_norm         residual normalization
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_pc_set_tolerance(cs_sles_pc_t  *pc,
                         double         precision,
                         double         r_norm)
{
  if (pc != NULL) {

    if (pc->context != NULL && pc->tolerance_func != NULL)
      pc->tolerance_func(pc->context, precision, r_norm);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Setup sparse linear equation preconditioner.
 *
 * Use of this function is optional: if a \ref cs_sles_solve is called
 * for the same system before this function is called, the latter will be
 * called automatically.
 *
 * If no options were previously provided for the matching system,
 * default options will be used.
 *
 * \param[in, out]  pc         pointer to preconditioner object
 * \param[in]       name       linear system name
 * \param[in]       a          matrix
 * \param[in]       accel      use accelerator version ?
 * \param[in]       verbosity  verbosity level
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_pc_setup(cs_sles_pc_t       *pc,
                 const char         *name,
                 const cs_matrix_t  *a,
                 bool                accel,
                 int                 verbosity)
{
  if (pc != NULL) {
    if (pc->context != NULL && pc->setup_func != NULL)
      pc->setup_func(pc->context, name, a, accel, verbosity);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply a preconditioner.
 *
 * If no options were previously provided for the matching system,
 * default options will be used.
 *
 * In cases where it is desired that the preconditioner modify a vector
 * "in place", x_in should be set to NULL, and x_out contain the vector to
 * be modified (\f$x_{out} \leftarrow M^{-1}x_{out})\f$).
 *
 * \param[in, out]  pc             pointer to preconditioner object
 * \param[in]       x_in           input vector
 * \param[in, out]  x_out          input/output vector
 *
 * \return  preconditioner application status
 */
/*----------------------------------------------------------------------------*/

cs_sles_pc_state_t
cs_sles_pc_apply(cs_sles_pc_t        *pc,
                 cs_real_t           *x_in,
                 cs_real_t           *x_out)
{
  return pc->apply_func(pc->context, x_in, x_out);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free preconditioner setup.
 *
 * This function frees resolution-related data, such as multigrid hierarchy,
 * preconditioning, and any other temporary arrays or objects required for
 * resolution, but maintains context information such as that used for
 * logging (especially performance data).
 *
 * \param[in, out]  pc  pointer to preconditioner object
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_pc_free(cs_sles_pc_t  *pc)
{
  if (pc != NULL) {
    if (pc->free_func != NULL)
      pc->free_func(pc->context);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log preconditioner setup, history and performance data.
 *
 * This function frees resolution-related data, such as multigrid hierarchy,
 * preconditioning, and any other temporary arrays or objects required for
 * resolution, but maintains context information such as that used for
 * logging (especially performance data).
 *
 * \param[in, out]  pc        pointer to preconditioner object
 * \param[in]       log_type  log type
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_pc_log(cs_sles_pc_t  *pc,
               cs_log_t       log_type)
{
  if (pc->log_func != NULL)
    pc->log_func(pc->context, log_type);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create an "identity" (or null) preconditioner.
 *
 * \return  pointer to newly created preconditioner object.
 */
/*----------------------------------------------------------------------------*/

cs_sles_pc_t *
cs_sles_pc_none_create(void)
{
  cs_sles_pc_poly_t *pcp = _sles_pc_poly_create();

  pcp->poly_degree = -1;

  cs_sles_pc_t *pc = cs_sles_pc_define(pcp,
                                       _sles_pc_poly_get_type,
                                       _sles_pc_poly_setup_none,
                                       NULL,
                                       _sles_pc_poly_apply_none,
                                       _sles_pc_poly_free,
                                       NULL,
                                       _sles_pc_poly_clone,
                                       _sles_pc_poly_destroy);

  return pc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a Jacobi preconditioner.
 *
 * \return  pointer to newly created preconditioner object.
 */
/*----------------------------------------------------------------------------*/

cs_sles_pc_t *
cs_sles_pc_jacobi_create(void)
{
  cs_sles_pc_poly_t *pcp = _sles_pc_poly_create();

  pcp->poly_degree = 0;

  cs_sles_pc_t *pc = cs_sles_pc_define(pcp,
                                       _sles_pc_poly_get_type,
                                       _sles_pc_poly_setup,
                                       NULL,
                                       _sles_pc_poly_apply_jacobi,
                                       _sles_pc_poly_free,
                                       NULL,
                                       _sles_pc_poly_clone,
                                       _sles_pc_poly_destroy);

  return pc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a Polynomial preconditioner of degree 1.
 *
 * \return  pointer to newly created preconditioner object.
 */
/*----------------------------------------------------------------------------*/

cs_sles_pc_t *
cs_sles_pc_poly_1_create(void)
{
  cs_sles_pc_poly_t *pcp = _sles_pc_poly_create();

  pcp->poly_degree = 1;

  cs_sles_pc_t *pc = cs_sles_pc_define(pcp,
                                       _sles_pc_poly_get_type,
                                       _sles_pc_poly_setup,
                                       NULL,
                                       _sles_pc_poly_apply_poly,
                                       _sles_pc_poly_free,
                                       NULL,
                                       _sles_pc_poly_clone,
                                       _sles_pc_poly_destroy);

  return pc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a Polynomial preconditioner of degree 2.
 *
 * \return  pointer to newly created preconditioner object.
 */
/*----------------------------------------------------------------------------*/

cs_sles_pc_t *
cs_sles_pc_poly_2_create(void)
{
  cs_sles_pc_poly_t *pcp = _sles_pc_poly_create();

  pcp->poly_degree = 2;

  cs_sles_pc_t *pc = cs_sles_pc_define(pcp,
                                       _sles_pc_poly_get_type,
                                       _sles_pc_poly_setup,
                                       NULL,
                                       _sles_pc_poly_apply_poly,
                                       _sles_pc_poly_free,
                                       NULL,
                                       _sles_pc_poly_clone,
                                       _sles_pc_poly_destroy);

  return pc;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
