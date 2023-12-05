/*============================================================================
 * Functions dedicated to the linear algebra settings and operations in case
 * of CDO cell-based schemes with a monolithic velocity-pressure coupling
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
#include "cs_cdocb_priv.h"
#include "cs_cdocb_scaleq.h"
#include "cs_fp_exception.h"
#include "cs_matrix_default.h"
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

#include "cs_cdocb_monolithic_sles.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
 * \file cs_cdocb_monolithic_sles.c
 *
 * \brief Functions dedicated to to the linear algebra settings and operations
 *        in case of CDO cell-based schemes with a monolithic coupling
 */

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_CDOCB_MONOLITHIC_SLES_DBG      0

/* Block size for superblock algorithm */

#define CS_SBLOCK_BLOCK_SIZE 60

/* Cache line multiple, in cs_real_t units */

#define CS_CL  (CS_CL_SIZE/8)

/* This structure follow notations given in the article entitled
 * "An iterative generalized Golub-Kahan algorithm for problems in structural
 *  mechanics" by M. Arioli, C. Kruse, U. Ruede and N. Tardieu
 *
 * M space is isomorphic to the face flux space (size = n_faces)
 * N space is isomorphic to the potential space (size = n_cells)
 */


/*============================================================================
 * Private variables
 *============================================================================*/

/* Pointer to shared structures */

static const cs_cdo_connect_t       *cs_shared_connect;
static const cs_cdo_quantities_t    *cs_shared_quant;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create an empty cs_cdocb_monolithic_sles_t structure
 *
 * \param[in] n_faces     number of faces (interior + border)
 * \param[in] n_cells     number of cells
 *
 * \return a pointer to a newly allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_cdocb_monolithic_sles_t *
cs_cdocb_monolithic_sles_create(cs_lnum_t    n_faces,
                                cs_lnum_t    n_cells)
{
  cs_cdocb_monolithic_sles_t  *msles = NULL;

  BFT_MALLOC(msles, 1, cs_cdocb_monolithic_sles_t);

  msles->div_op = NULL;

  msles->graddiv_coef = 0.;

  msles->sles = NULL;
  msles->schur_sles = NULL;

  msles->n_faces = n_faces;
  msles->n_cells = n_cells;

  msles->flux = NULL;
  msles->potential = NULL;

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
cs_cdocb_monolithic_sles_clean(cs_cdocb_monolithic_sles_t   *msles)
{
  if (msles == NULL)
    return;

  cs_sles_free(msles->sles);
  cs_sles_free(msles->schur_sles);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free memory related to cs_cdocb_monolithic_sles_t structure
 *
 * \param[in, out]  p_msles  double pointer to the structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_monolithic_sles_free(cs_cdocb_monolithic_sles_t   **p_msles)
{
  cs_cdocb_monolithic_sles_t  *msles = *p_msles;

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
cs_cdocb_monolithic_sles_init_sharing(const cs_cdo_connect_t        *connect,
                                      const cs_cdo_quantities_t     *quant)
{
  /* Assign static const pointers */

  cs_shared_connect = connect;
  cs_shared_quant = quant;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free if needed structure(s) associated CDO cell-based schemes with
 *         a monolithic velocity-pressure coupling
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_monolithic_sles_finalize(void)
{

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Start setting-up the solver when a monolithic
 *         algorithm is used to couple the system.
 *         No mesh information is available at this stage.
 *
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_monolithic_set_sles(cs_equation_param_t   *eqp,
                             void                  *context)
{
  cs_param_sles_t  *slesp = eqp->sles_param;
  cs_cdocb_scaleq_t *eqc = (cs_cdocb_scaleq_t *) context;

  int  field_id = eqc->var_field_id;

  slesp->field_id = field_id;

#if defined(HAVE_MUMPS)
  if (slesp->solver == CS_PARAM_ITSOL_MUMPS)
    cs_sles_mumps_define(field_id,
                         NULL,
                         slesp,
                         cs_user_sles_mumps_hook,
                         NULL);
#endif  /* HAVE_MUMPS */

  /* Define the level of verbosity for SLES structure */

  if (slesp->verbosity > 1) {

    cs_sles_t  *sles = cs_sles_find_or_add(field_id, NULL);

    /* Set verbosity */

    cs_sles_set_verbosity(sles, slesp->verbosity);

    eqc->msles->sles = sles;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve a linear system arising from the discretization of the
 *         diffsuion equation with a CDO cell-based approach.
 *         The full system is treated as one block and solved as it is.
 *         In this situation, PETSc or MUMPS are usually considered.
 *
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in]      sh       pointer to a cs_cdo_system_helper_t structure
 * \param[in, out] msles    pointer to a cs_cdocb_monolithic_sles_t structure
 *
 * \return the (cumulated) number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdocb_monolithic_solve(const cs_equation_param_t     *eqp,
                          const cs_cdo_system_helper_t  *sh,
                          cs_cdocb_monolithic_sles_t    *msles)
{
  assert(sh != NULL);
  assert(sh->n_blocks == 1);

  const cs_range_set_t  *range_set = cs_cdo_system_get_range_set(sh, 0);
  const cs_matrix_t  *matrix = cs_cdo_system_get_matrix(sh, 0);
  const cs_lnum_t  n_cols = cs_matrix_get_n_columns(matrix);
  const cs_lnum_t  n_faces = msles->n_faces;
  const cs_lnum_t  n_cells = msles->n_cells;
  const cs_lnum_t  n_scatter_elts = n_faces + n_cells;

  /* De-interlace the velocity array and the rhs for the face DoFs */

  cs_real_t  *sol = NULL;
  BFT_MALLOC(sol, CS_MAX(n_cols, n_scatter_elts), cs_real_t);

  cs_real_t  *b = NULL;
  BFT_MALLOC(b, n_scatter_elts, cs_real_t);

# pragma omp parallel for if (CS_THR_MIN > n_faces)     \
  shared(msles, sol, b) firstprivate(n_faces)
  for (cs_lnum_t f = 0; f < n_faces; f++) {

    sol[f] = msles->flux[f];

    b[f] = sh->rhs[f];

  }

  /* Add the pressure related elements */

  memcpy(sol + n_faces, msles->potential, n_cells*sizeof(cs_real_t));
  memcpy(b + n_faces, sh->rhs + n_faces, n_cells*sizeof(cs_real_t));

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

  const cs_param_sles_t  *sles_param = eqp->sles_param;
  const double  r_norm = 1.0; /* No renormalization by default (TODO) */

  cs_real_t  rtol = sles_param->cvg_param.rtol;
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

  if (sles_param->verbosity > 1)
    cs_log_printf(CS_LOG_DEFAULT, "  <%20s/sles_cvg_code=%-d> n_iters %d |"
                  " residual % -8.4e\n",
                  eqp->name, code, n_iters, residual);

  /* sol is computed and stored in a "gather" view. Switch to a "scatter"
     view */

  cs_range_set_scatter(range_set,
                       CS_REAL_TYPE, 1, /* type and stride */
                       sol, sol);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOCB_MONOLITHIC_SLES_DBG > 1
  cs_range_set_scatter(range_set,
                       CS_REAL_TYPE, 1, /* type and stride */
                       b, b);

  cs_dbg_fprintf_system(eqp->name,
                        -1,
                        CS_CDOCB_MONOLITHIC_SLES_DBG,
                        sol, b, n_faces);
#endif

  /* Switch from sol (not interlaced) to flux and potential */

  cs_real_t  *flux = msles->flux;

  /* Copy the part of the solution array related to the pressure in cells */

  memcpy(msles->potential, sol + n_faces, n_cells*sizeof(cs_real_t));

# pragma omp parallel for if (CS_THR_MIN > n_faces)     \
  shared(msles, sol) firstprivate(n_faces)
  for (cs_lnum_t f = 0; f < n_faces; f++)
    flux[f] = sol[f];

  /* Free what can be freed at this stage */

  BFT_FREE(sol);
  BFT_FREE(b);

  return n_iters;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
