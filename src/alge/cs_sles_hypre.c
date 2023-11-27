/*============================================================================
 * Sparse Linear Equation Solvers using HYPRE
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
 * HYPRE headers
 *----------------------------------------------------------------------------*/

/* Avoid warnings due to previous values */
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_URL
#undef PACKAGE_VERSION

#include <HYPRE_config.h>
#include <HYPRE_krylov.h>
#include <HYPRE_parcsr_ls.h>
#include <HYPRE_utilities.h>

#if !defined(HYPRE_RELEASE_NUMBER)
#define HYPRE_RELEASE_NUMBER 0
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_base_accel.h"
#include "cs_log.h"
#include "cs_fp_exception.h"
#include "cs_halo.h"
#include "cs_matrix.h"
#include "cs_matrix_default.h"
#include "cs_matrix_hypre.h"
#include "cs_matrix_hypre_priv.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_sles.h"
#include "cs_sles_hypre.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_sles_hypre.c

  \brief handling of HYPRE-based linear solvers

  \page sles_hypre HYPRE-based linear solvers.

  \typedef cs_sles_hypre_setup_hook_t

  \brief Function pointer for settings of a HYPRE solver setup.

  This function is called during the setup stage for a HYPRE solver.

  When first called, the solver argument is NULL, and must be created
  using HYPRE functions.

  Note: if the context pointer is non-NULL, it must point to valid data
  when the selection function is called so that value or structure should
  not be temporary (i.e. local);

  \param[in]       verbosity  verbosity level
  \param[in, out]  context    pointer to optional (untyped) value or structure
  \param[in, out]  solver     handle to HYPRE solver (to be cast as HYPRE_Solver)
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

typedef struct _cs_sles_hypre_setup_t {

  cs_matrix_coeffs_hypre_t *coeffs;  /* HYPRE matrix and vectors */

  HYPRE_Solver          solver;            /* Solver data */
  HYPRE_Solver          precond;           /* Preconditioner data */

} cs_sles_hypre_setup_t;

typedef struct _cs_sles_hypre_t {

  cs_sles_hypre_type_t  solver_type;       /* Solver type */
  cs_sles_hypre_type_t  precond_type;      /* Preconditioner type */

  int                   use_device;        /* O for host, 1 for device */

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
  cs_sles_hypre_setup_hook_t  *setup_hook;     /* Post setup function */

  cs_sles_hypre_setup_t     *setup_data;

} cs_sles_hypre_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

static int  _n_hypre_systems = 0;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return name of hypre solver type.
 *
 * \param[in]      solver_type   solver type id
 *
 * \return  name od associated solver type.
 */
/*----------------------------------------------------------------------------*/

static const char *
_cs_hypre_type_name(cs_sles_hypre_type_t  solver_type)
{
  switch(solver_type) {

  case CS_SLES_HYPRE_BOOMERAMG:
    return N_("BoomerAMG");
    break;
  case CS_SLES_HYPRE_HYBRID:
    return N_("Hybrid");
    break;
  case CS_SLES_HYPRE_ILU:
    return N_("ILU");
    break;
  case CS_SLES_HYPRE_BICGSTAB:
    return N_("BiCGSTAB");
    break;
  case CS_SLES_HYPRE_GMRES:
    return N_("GMRES");
    break;
  case CS_SLES_HYPRE_FLEXGMRES:
    return N_("Flexible GMRES");
    break;
  case CS_SLES_HYPRE_LGMRES:
    return N_("LGMRES");
    break;
  case CS_SLES_HYPRE_PCG:
    return N_("PCG");
    break;
  case CS_SLES_HYPRE_EUCLID:
    return N_("EUCLID");
    break;
  case CS_SLES_HYPRE_PARASAILS:
    return N_("ParaSails");
    break;
  case CS_SLES_HYPRE_NONE:
    return N_("None");
    break;

  default:
    return NULL;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Ensure MPI is initialized if present and needed
 */
/*----------------------------------------------------------------------------*/

static void
_ensure_mpi_init(void)
{
#if defined(HAVE_MPI) && defined(HYPRE_HAVE_MPI)
  int flag = 0;
  MPI_Initialized(&flag);
  if (!flag) {
    int mpi_threads;
    MPI_Init_thread(NULL, NULL, MPI_THREAD_FUNNELED, &mpi_threads);
  }
#endif
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define and associate a HYPRE linear system solver
 *        for a given field or equation name.
 *
 * If this system did not previously exist, it is added to the list of
 * "known" systems. Otherwise, its definition is replaced by the one
 * defined here.
 *
 * This is a utility function: if finer control is needed, see
 * \ref cs_sles_define and \ref cs_sles_petsc_create.
 *
 * The associated solver required that the matrix passed to it is a HYPRE
 * matrix (see cs_matrix_set_type_hypre).
 *
 * Note that this function returns a pointer directly to the iterative solver
 * management structure. This may be used to set further options.
 * If needed, \ref cs_sles_find may be used to obtain a pointer to the matching
 * \ref cs_sles_t container.
 *
 * \param[in]      f_id          associated field id, or < 0
 * \param[in]      name          associated name if f_id < 0, or NULL
 * \param[in]      solver_type   HYPRE solver type
 * \param[in]      precond_type  HYPRE preconditioner type
 * \param[in]      setup_hook    pointer to optional setup epilogue function
 * \param[in,out]  context       pointer to optional (untyped) value or
 *                               structure for setup_hook, or NULL
 *
 * \return  pointer to newly created iterative solver info object.
 */
/*----------------------------------------------------------------------------*/

cs_sles_hypre_t *
cs_sles_hypre_define(int                          f_id,
                     const char                  *name,
                     cs_sles_hypre_type_t         solver_type,
                     cs_sles_hypre_type_t         precond_type,
                     cs_sles_hypre_setup_hook_t  *setup_hook,
                     void                        *context)
{
  if (solver_type < 0 || solver_type >= CS_SLES_HYPRE_NONE)
    bft_error(__FILE__, __LINE__, 0,
              _("Incorrect solver type argument %d for HYPRE."),
              (int)solver_type);

  cs_sles_hypre_t *c = cs_sles_hypre_create(solver_type,
                                            precond_type,
                                            setup_hook,
                                            context);

  cs_sles_t *sc = cs_sles_define(f_id,
                                 name,
                                 c,
                                 "cs_sles_hypre_t",
                                 cs_sles_hypre_setup,
                                 cs_sles_hypre_solve,
                                 cs_sles_hypre_free,
                                 cs_sles_hypre_log,
                                 cs_sles_hypre_copy,
                                 cs_sles_hypre_destroy);

  cs_sles_set_error_handler(sc,
                            cs_sles_hypre_error_post_and_abort);

  return c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create HYPRE linear system solver info and context.
 *
 * In case of rotational periodicity for a block (non-scalar) matrix,
 * the matrix type will be forced to MATSHELL ("shell") regardless
 * of the option used.
 *
 * \param[in]  solver_type   HYPRE solver type
 * \param[in]  precond_type  HYPRE preconditioner type
 * \param[in]  setup_hook    pointer to optional setup epilogue function
 * \param[in]  context       pointer to optional (untyped) value or structure
 *                           for setup_hook, or NULL
 *
 * \return  pointer to newly created linear system object.
 */
/*----------------------------------------------------------------------------*/

cs_sles_hypre_t *
cs_sles_hypre_create(cs_sles_hypre_type_t         solver_type,
                     cs_sles_hypre_type_t         precond_type,
                     cs_sles_hypre_setup_hook_t  *setup_hook,
                     void                        *context)
{
  cs_sles_hypre_t *c;

  if (_n_hypre_systems == 0) {
    _ensure_mpi_init();
#if HYPRE_RELEASE_NUMBER < 22900
    HYPRE_Init();
#else
    if (HYPRE_Initialized() == 0)
      HYPRE_Initialize();
#endif
  }
  _n_hypre_systems += 1;


  BFT_MALLOC(c, 1, cs_sles_hypre_t);

  c->solver_type = solver_type;
  c->precond_type = precond_type;

  int device_id = cs_get_device_id();
  c->use_device = (device_id < 0) ? 0 : 1;

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
  c->setup_data = NULL;

  return c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy iterative sparse linear system solver info and context.
 *
 * \param[in, out]  context  pointer to iterative solver info and context
 *                           (actual type: cs_sles_hypre_t  **)
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_hypre_destroy(void **context)
{
  cs_sles_hypre_t *c = (cs_sles_hypre_t *)(*context);
  if (c != NULL) {
    cs_sles_hypre_free(c);
    BFT_FREE(c);
    *context = c;
  }

  _n_hypre_systems -= 1;
  if (_n_hypre_systems == 0) {
    HYPRE_Finalize();
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create HYPRE sparse linear system solver info and context
 *        based on existing info and context.
 *
 * \param[in]  context  pointer to reference info and context
 *                     (actual type: cs_sles_hypre_t  *)
 *
 * \return  pointer to newly created solver info object.
 *          (actual type: cs_sles_hypre_t  *)
 */
/*----------------------------------------------------------------------------*/

void *
cs_sles_hypre_copy(const void  *context)
{
  cs_sles_hypre_t *d = NULL;

  if (context != NULL) {
    const cs_sles_hypre_t *c = context;
    d = cs_sles_hypre_create(c->solver_type,
                             c->precond_type,
                             c->setup_hook,
                             c->hook_context);
  }

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Error handler for HYPRE solver.
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
cs_sles_hypre_error_post_and_abort(cs_sles_t                    *sles,
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

  const char *name = cs_sles_get_name(sles);

  const cs_sles_hypre_t  *c = cs_sles_get_context(sles);
  CS_UNUSED(c);

  const char *error_type[] = {N_("divergence"), N_("breakdown")};
  int err_id = (state == CS_SLES_BREAKDOWN) ? 1 : 0;

  bft_error(__FILE__, __LINE__, 0,
            _("HYPRE: error (%s) solving for %s"),
            _(error_type[err_id]),
            name);

  return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log sparse linear equation solver info.
 *
 * \param[in]  context   pointer to iterative solver info and context
 *                       (actual type: cs_sles_hypre_t  *)
 * \param[in]  log_type  log type
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_hypre_log(const void  *context,
                  cs_log_t     log_type)
{
  const cs_sles_hypre_t  *c = context;

  if (log_type == CS_LOG_SETUP) {

    cs_log_printf(log_type,
                  _("  Solver type:                       HYPRE (%s)\n"),
                  _cs_hypre_type_name(c->solver_type));
    if (c->precond_type < CS_SLES_HYPRE_NONE)
      cs_log_printf(log_type,
                    _("    Preconditioning:                 %s\n"),
                    _cs_hypre_type_name(c->precond_type));
    if (c->use_device)
      cs_log_printf(log_type,
                    _("    Accelerated device:              enabled\n"));

  }

  else if (log_type == CS_LOG_PERFORMANCE) {

    int n_calls = c->n_solves;
    int n_it_min = c->n_iterations_min;
    int n_it_max = c->n_iterations_max;
    int n_it_mean = 0;

    if (n_calls > 0)
      n_it_mean = (int)(  c->n_iterations_tot
          / ((int long long)n_calls));

    cs_log_printf(log_type,
                  _("\n"
                    "  Solver type:                   HYPRE (%s)\n"),
                  _cs_hypre_type_name(c->solver_type));
    if (c->precond_type < CS_SLES_HYPRE_NONE)
      cs_log_printf(log_type,
                    _("    Preconditioning:             %s\n"),
                    _cs_hypre_type_name(c->precond_type));
    if (c->use_device)
      cs_log_printf(log_type,
                    _("    Accelerated device:          enabled\n"));
    cs_log_printf(log_type,
                  _("  Number of setups:              %12d\n"
                    "  Number of calls:               %12d\n"
                    "  Minimum number of iterations:  %12d\n"
                    "  Maximum number of iterations:  %12d\n"
                    "  Mean number of iterations:     %12d\n"
                    "  Construction:                  %12.3f\n"
                    "  Resolution:                    %12.3f\n"),
                  c->n_setups, n_calls, n_it_min, n_it_max, n_it_mean,
                  c->t_setup.nsec*1e-9,
                  c->t_solve.nsec*1e-9);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Setup iterative sparse linear equation solver.
 *
 * \param[in, out]  context    pointer to iterative solver info and context
 *                             (actual type: cs_sles_hypre_t  *)
 * \param[in]       name       pointer to system name
 * \param[in]       a          associated matrix
 * \param[in]       verbosity  associated verbosity
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_hypre_setup(void               *context,
                    const char         *name,
                    const cs_matrix_t  *a,
                    int                 verbosity)
{
  CS_UNUSED(name);

  cs_timer_t t0;
  t0 = cs_timer_time();

  cs_sles_hypre_t  *c  = context;
  cs_sles_hypre_setup_t *sd = c->setup_data;
  if (sd == NULL) {
    BFT_MALLOC(c->setup_data, 1, cs_sles_hypre_setup_t);
    sd = c->setup_data;
    sd->solver = NULL;
    sd->precond = NULL;
  }

  const char expected_matrix_type[] = "HYPRE_PARCSR";
  if (strcmp(cs_matrix_get_type_name(a), expected_matrix_type))
    bft_error(__FILE__, __LINE__, 0,
              _("HYPRE [%s]:\n"
                "  expected matrix type: %s\n"
                "  provided matrix type: %s"),
              name, expected_matrix_type, cs_matrix_get_type_name(a));

  sd->coeffs = cs_matrix_hypre_get_coeffs(a);

  MPI_Comm comm = cs_glob_mpi_comm;
  if (comm == MPI_COMM_NULL)
    comm = MPI_COMM_WORLD;

  bool have_set_pc = true;
  HYPRE_PtrToParSolverFcn solve_ftn[2] = {NULL, NULL};
  HYPRE_PtrToParSolverFcn setup_ftn[2] = {NULL, NULL};

  for (int i = 0; i < 2; i++) {

    HYPRE_Solver hs = (i == 0) ? sd->precond : sd->solver;
    cs_sles_hypre_type_t hs_type = (i == 0) ? c->precond_type : c->solver_type;

    /* hs should be NULL at this point, unless we do not relly free
     it (when cs_sles_hypre_free is called (for example to amortize setup) */

    if (hs != NULL || hs_type >= CS_SLES_HYPRE_NONE)
      continue;

    switch(hs_type) {

    case CS_SLES_HYPRE_BOOMERAMG:
      {
        HYPRE_BoomerAMGCreate(&hs);

        if (verbosity > 2 ) {
          HYPRE_BoomerAMGSetPrintLevel(hs, 1);
          HYPRE_BoomerAMGSetPrintFileName(hs, "hypre.log");
        }

        /* Default settings for device */
        if (c->use_device == 1) {
          HYPRE_BoomerAMGSetRelaxType(hs, 18);      /* 3, 4, 6, 7, 18, 11, 12 */
          HYPRE_BoomerAMGSetRelaxOrder(hs, 0);      /* must be false */
          HYPRE_BoomerAMGSetCoarsenType(hs, 8);     /* PMIS */
          HYPRE_BoomerAMGSetInterpType(hs, 18);     /* 3, 15, 6, 14, 17, 18 */
          HYPRE_BoomerAMGSetAggInterpType(hs, 7);   /* 5 or 7 */
          HYPRE_BoomerAMGSetKeepTranspose(hs, 1);   /* keep transpose to
                                                       avoid SpMTV */
          HYPRE_BoomerAMGSetRAP2(hs, 0);            /* RAP in two multiplications
                                                       (default: FALSE) */

          HYPRE_BoomerAMGSetAggNumLevels(hs, 2);
          HYPRE_BoomerAMGSetPMaxElmts(hs, 4);

          //HYPRE_SetSpGemmUseVendor(0);
          HYPRE_BoomerAMGSetAggNumLevels(hs, 2);
        }

        /* Default settings for host */
        else if (c->use_device == 0) {

          HYPRE_BoomerAMGSetCoarsenType(hs, 10); /* HMIS */
          HYPRE_BoomerAMGSetPMaxElmts(hs, 4);
          HYPRE_BoomerAMGSetInterpType(hs, 7);   /* extended+e */
          HYPRE_BoomerAMGSetRelaxType(hs, 6);    /* Sym G.S./Jacobi hybrid */
          HYPRE_BoomerAMGSetRelaxOrder(hs, 0);

        }

        /* Defaults for both host and device */
        HYPRE_BoomerAMGSetAggNumLevels(hs, 2);
        HYPRE_BoomerAMGSetStrongThreshold(hs, 0.5);   /* 2d=>0.25 3d=>0.5 */

        solve_ftn[i] = HYPRE_BoomerAMGSolve;
        setup_ftn[i] = HYPRE_BoomerAMGSetup;

        if (i == 0) { /* preconditioner */
          HYPRE_BoomerAMGSetTol(hs, 0.0);
          HYPRE_BoomerAMGSetMaxIter(hs, 1);
        }
        else { /* solver */
          have_set_pc = false;
          HYPRE_BoomerAMGSetMaxIter(hs, 1000);
        }

      }
      break;

    case CS_SLES_HYPRE_HYBRID:
      {
        HYPRE_ParCSRHybridCreate(&hs);

        if (verbosity > 2 ) {
          HYPRE_ParCSRHybridSetPrintLevel(hs, 2);  /* Print solve info */
          HYPRE_ParCSRHybridSetLogging(hs, 1);     /* Needed to get info later */
        }

        solve_ftn[i] = HYPRE_ParCSRHybridSolve;
        setup_ftn[i] = HYPRE_ParCSRHybridSetup;

        if (i == 1 && solve_ftn[0] != NULL) {  /* solver */
          HYPRE_ParCSRHybridSetPrecond(hs,
                                       solve_ftn[0],
                                       setup_ftn[0],
                                       sd->precond);
        }
      }
      break;

    case CS_SLES_HYPRE_ILU:
      {
        HYPRE_ILUCreate(&hs);

        solve_ftn[i] = HYPRE_ILUSolve;
        setup_ftn[i] = HYPRE_ILUSetup;

        if (i == 0) { /* preconditioner */
          HYPRE_ILUSetTol(hs, 0.);
        }
        else { /* solver */
          have_set_pc = false;
        }

        if (verbosity > 2 ) {
          HYPRE_ILUSetPrintLevel(hs, 2);  /* Print solve info */
          HYPRE_ILUSetLogging(hs, 1);     /* Needed to get info later */
        }
      }
      break;

    case CS_SLES_HYPRE_BICGSTAB:
      {
        HYPRE_ParCSRBiCGSTABCreate(comm, &hs);

        if (verbosity > 2 ) {
          HYPRE_BiCGSTABSetPrintLevel(hs, 2);  /* Print solve info */
          HYPRE_BiCGSTABSetLogging(hs, 1);     /* Needed to get run info later */
        }

        solve_ftn[i] = HYPRE_ParCSRBiCGSTABSolve;
        setup_ftn[i] = HYPRE_ParCSRBiCGSTABSetup;

        if (i == 0) { /* preconditioner */
          HYPRE_ParCSRBiCGSTABSetTol(hs, 0.0);
          HYPRE_ParCSRBiCGSTABSetMaxIter(hs, 1);
          solve_ftn[i] = HYPRE_ParCSRBiCGSTABSolve;
          setup_ftn[i] = HYPRE_ParCSRBiCGSTABSetup;
        }
        else { /* solver */
          HYPRE_BiCGSTABSetMaxIter(hs, 1000);  /* Max iterations */
          if (solve_ftn[0] != NULL)
            HYPRE_ParCSRBiCGSTABSetPrecond(hs,
                                           solve_ftn[0],
                                           setup_ftn[0],
                                           sd->precond);
        }
      }
      break;

    case CS_SLES_HYPRE_GMRES:
      {
        HYPRE_ParCSRGMRESCreate(comm, &hs);

        if (verbosity > 2 ) {
          HYPRE_GMRESSetPrintLevel(hs, 2);  /* Print solve info */
          HYPRE_GMRESSetLogging(hs, 1);     /* Needed to get run info later */
          solve_ftn[i] = HYPRE_ParCSRGMRESSolve;
          setup_ftn[i] = HYPRE_ParCSRGMRESSetup;
        }

        solve_ftn[i] = HYPRE_ParCSRGMRESSolve;
        setup_ftn[i] = HYPRE_ParCSRGMRESSetup;

        if (i == 0) { /* preconditioner */
          HYPRE_ParCSRGMRESSetTol(hs, 0.0);
          HYPRE_ParCSRGMRESSetMaxIter(hs, 1);
        }
        else { /* solver */
          HYPRE_GMRESSetMaxIter(hs, 1000);  /* Max iterations */
          if (solve_ftn[0] != NULL)
            HYPRE_ParCSRGMRESSetPrecond(hs,
                                        solve_ftn[0],
                                        setup_ftn[0],
                                        sd->precond);
        }
      }
      break;

    case CS_SLES_HYPRE_FLEXGMRES:
      {
        HYPRE_ParCSRFlexGMRESCreate(comm, &hs);

        if (verbosity > 2 ) {
          HYPRE_FlexGMRESSetPrintLevel(hs, 2);  /* Print solve info */
          HYPRE_FlexGMRESSetLogging(hs, 1);     /* Needed to get run info later */
        }

        solve_ftn[i] = HYPRE_ParCSRFlexGMRESSolve;
        setup_ftn[i] = HYPRE_ParCSRFlexGMRESSetup;

        if (i == 0) { /* preconditioner */
          HYPRE_ParCSRFlexGMRESSetTol(hs, 0.0);
          HYPRE_ParCSRFlexGMRESSetMaxIter(hs, 1);
        }
        else { /* solver */
          HYPRE_FlexGMRESSetMaxIter(hs, 1000);  /* Max iterations */
          if (solve_ftn[0] != NULL)
            HYPRE_ParCSRFlexGMRESSetPrecond(hs,
                                            solve_ftn[0],
                                            setup_ftn[0],
                                            sd->precond);
        }
      }
      break;

    case CS_SLES_HYPRE_LGMRES:
      {
        HYPRE_ParCSRLGMRESCreate(comm, &hs);

        if (verbosity > 2 ) {
          HYPRE_LGMRESSetPrintLevel(hs, 2);  /* Print solve info */
          HYPRE_LGMRESSetLogging(hs, 1);     /* Needed to get run info later */
        }

        solve_ftn[i] = HYPRE_ParCSRLGMRESSolve;
        setup_ftn[i] = HYPRE_ParCSRLGMRESSetup;

        if (i == 0) { /* preconditioner */
          HYPRE_ParCSRLGMRESSetTol(hs, 0.0);
          HYPRE_ParCSRLGMRESSetMaxIter(hs, 1);
        }
        else { /* solver */
          HYPRE_LGMRESSetMaxIter(hs, 1000);  /* Max iterations */
          if (solve_ftn[0] != NULL)
            HYPRE_ParCSRLGMRESSetPrecond(hs,
                                         solve_ftn[0],
                                         setup_ftn[0],
                                         sd->precond);
        }
      }
      break;

    case CS_SLES_HYPRE_PCG:
      {
        HYPRE_ParCSRPCGCreate(comm, &hs);

        if (verbosity > 2 ) {
          HYPRE_PCGSetPrintLevel(hs, 2);  /* Print solve info */
          HYPRE_PCGSetLogging(hs, 1);     /* Needed to get run info later */
        }

        solve_ftn[i] = HYPRE_ParCSRPCGSolve;
        setup_ftn[i] = HYPRE_ParCSRPCGSetup;

        if (i == 0) { /* preconditioner */
          HYPRE_ParCSRPCGSetTol(hs, 0.0);
          HYPRE_ParCSRPCGSetMaxIter(hs, 1);
        }
        else { /* solver */
          HYPRE_PCGSetMaxIter(hs, 1000);  /* Max iterations */
          if (solve_ftn[0] != NULL)
            HYPRE_ParCSRPCGSetPrecond(hs,
                                      solve_ftn[0],
                                      setup_ftn[0],
                                      sd->precond);
        }
      }
      break;

    case CS_SLES_HYPRE_EUCLID:
      {
        HYPRE_EuclidCreate(comm, &hs);

        solve_ftn[i] = HYPRE_EuclidSolve;
        setup_ftn[i] = HYPRE_EuclidSetup;

        if (i > 0) /* solver */
          bft_error(__FILE__, __LINE__, 0,
                    _("HYPRE: type (%s) is a preconditioner, not a solver."),
                    _cs_hypre_type_name(c->solver_type));

      }
      break;

    case CS_SLES_HYPRE_PARASAILS:
      {
        HYPRE_ParCSRParaSailsCreate(comm, &hs);

        if (verbosity > 2 ) {
          HYPRE_ParCSRParaSailsSetLogging(hs, 1);     /* Needed to get run
                                                         info later */
        }

        solve_ftn[i] = HYPRE_ParCSRParaSailsSolve;
        setup_ftn[i] = HYPRE_ParCSRParaSailsSetup;
        if (i > 0) /* solver */
          bft_error(__FILE__, __LINE__, 0,
                    _("HYPRE: type (%s) is a preconditioner, not a solver."),
                    _cs_hypre_type_name(c->solver_type));
      }
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                _("HYPRE: solver type (%s) not currently handled."),
                _cs_hypre_type_name(c->solver_type));

    }

    if (i == 0)
      sd->precond = hs;
    else
      sd->solver = hs;
  }

  if (sd->precond != NULL && have_set_pc == false)
    bft_error(__FILE__, __LINE__, 0,
              _("HYPRE: solver (%s) will ignore preconditioner (%s)."),
              _cs_hypre_type_name(c->solver_type),
              _cs_hypre_type_name(c->precond_type));

  /* Call optional setup hook for user setting changes */

  if (c->setup_hook != NULL)
    c->setup_hook(verbosity, c->hook_context, sd->solver);

  /* Now setup systems (where rhs and vx values may be different
     when solving, but their shapes and addresses are the same) */

  HYPRE_ParCSRMatrix par_a;              /* Associated matrix */
  HYPRE_ParVector p_x, p_rhs;

  HYPRE_IJMatrixGetObject(sd->coeffs->hm, (void **)&par_a);
  HYPRE_IJVectorGetObject(sd->coeffs->hx, (void **)&p_x);
  HYPRE_IJVectorGetObject(sd->coeffs->hy, (void **)&p_rhs);

  if (setup_ftn[1] != NULL)
    setup_ftn[1](sd->solver, par_a, p_rhs, p_x);
  else
    bft_error(__FILE__, __LINE__, 0,
              _("HYPRE: setup function for solver type (%s) not set."),
              _cs_hypre_type_name(c->solver_type));

  /* Update return values */

  c->n_setups += 1;

  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(c->t_setup), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Call HYPRE linear equation solver.
 *
 * \param[in, out]  context        pointer to iterative solver info and context
 *                                 (actual type: cs_sles_hypre_t  *)
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
cs_sles_hypre_solve(void                *context,
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

  cs_timer_t t0;
  t0 = cs_timer_time();
  cs_sles_convergence_state_t cvg = CS_SLES_ITERATING;
  cs_sles_hypre_t  *c = context;
  cs_sles_hypre_setup_t  *sd = c->setup_data;
  cs_lnum_t n_rows = cs_matrix_get_n_rows(a);
  double res;
  HYPRE_Int its, max_its;
  HYPRE_Int ierr = -1;

  if (sd == NULL) {
    cs_sles_hypre_setup(c, name, a, verbosity);
    sd = c->setup_data;
  }

  HYPRE_ParCSRMatrix par_a;              /* Associted matrix */
  HYPRE_IJMatrixGetObject(sd->coeffs->hm, (void **)&par_a);

  cs_alloc_mode_t  amode = CS_ALLOC_HOST;
  if (sd->coeffs->memory_location != HYPRE_MEMORY_HOST)
    amode = CS_ALLOC_HOST_DEVICE_SHARED;

  HYPRE_Real *_t = NULL;

  /* Set RHS and starting solution */

  if (sizeof(cs_real_t) == sizeof(HYPRE_Real) && amode == CS_ALLOC_HOST) {
    HYPRE_IJVectorSetValues(sd->coeffs->hx, n_rows, NULL, vx);
    HYPRE_IJVectorSetValues(sd->coeffs->hy, n_rows, NULL, rhs);
  }
  else {
    CS_MALLOC_HD(_t, n_rows, HYPRE_Real, amode);
    for (HYPRE_BigInt ii = 0; ii < n_rows; ii++) {
      _t[ii] = vx[ii];
    }
    HYPRE_IJVectorSetValues(sd->coeffs->hx, n_rows, NULL, _t);
    for (HYPRE_BigInt ii = 0; ii < n_rows; ii++) {
      _t[ii] = rhs[ii];
    }
    HYPRE_IJVectorSetValues(sd->coeffs->hy, n_rows, NULL, _t);
  }

  HYPRE_IJVectorAssemble(sd->coeffs->hx);
  HYPRE_IJVectorAssemble(sd->coeffs->hy);

  HYPRE_ParVector p_x, p_rhs;
  HYPRE_IJVectorGetObject(sd->coeffs->hx, (void **)&p_x);
  HYPRE_IJVectorGetObject(sd->coeffs->hy, (void **)&p_rhs);

  switch(c->solver_type) {

  case CS_SLES_HYPRE_BOOMERAMG:
    {
      /* Finalize setup and solve; no absolute tolerance is available,
         so we use the available function. */
      HYPRE_BoomerAMGSetTol(sd->solver, precision*r_norm);
      ierr = HYPRE_BoomerAMGSolve(sd->solver, par_a, p_rhs, p_x);

      /* Get solution and information */
      HYPRE_BoomerAMGGetFinalRelativeResidualNorm(sd->solver, &res);
      HYPRE_BoomerAMGGetNumIterations(sd->solver, &its);
    }
    break;

  case CS_SLES_HYPRE_HYBRID:
    {
      /* Finalize setup and solve */
      HYPRE_ParCSRHybridSetAbsoluteTol(sd->solver, precision*r_norm);
      ierr = HYPRE_ParCSRHybridSolve(sd->solver, par_a, p_rhs, p_x);

      /* Get solution and information */
      HYPRE_ParCSRHybridGetFinalRelativeResidualNorm(sd->solver, &res);
      HYPRE_ParCSRHybridGetNumIterations(sd->solver, &its);
    }
    break;

  case CS_SLES_HYPRE_ILU:
    {
      /* Finalize setup and solve */
      HYPRE_ILUSetTol(sd->solver, precision*r_norm);
      ierr = HYPRE_ILUSolve(sd->solver, par_a, p_rhs, p_x);

      /* Get solution and information */
      HYPRE_ILUGetFinalRelativeResidualNorm(sd->solver, &res);
      HYPRE_ILUGetNumIterations(sd->solver, &its);

      /* Manage feedback */
      if (ierr != 0 && HYPRE_CheckError(ierr, HYPRE_ERROR_CONV))
        cvg = CS_SLES_DIVERGED;
    }
    break;

  case CS_SLES_HYPRE_BICGSTAB:
    {
      /* Finalize setup and solve */
      HYPRE_BiCGSTABSetAbsoluteTol(sd->solver, precision*r_norm);
      HYPRE_BiCGSTABSetTol(sd->solver, 0);
      ierr = HYPRE_ParCSRBiCGSTABSolve(sd->solver, par_a, p_rhs, p_x);

      /* Get solution and information */
      HYPRE_ParCSRBiCGSTABGetFinalRelativeResidualNorm(sd->solver, &res);
      HYPRE_ParCSRBiCGSTABGetNumIterations(sd->solver, &its);
    }
    break;

  case CS_SLES_HYPRE_GMRES:
    {
      /* Finalize setup and solve */
      HYPRE_GMRESSetAbsoluteTol(sd->solver, precision*r_norm);
      HYPRE_GMRESSetTol(sd->solver, 0);
      ierr = HYPRE_ParCSRGMRESSolve(sd->solver, par_a, p_rhs, p_x);

      /* Get solution and information */
      HYPRE_ParCSRGMRESGetFinalRelativeResidualNorm(sd->solver, &res);
      HYPRE_ParCSRGMRESGetNumIterations(sd->solver, &its);

      /* Manage feedback */
      if (ierr != 0 && HYPRE_CheckError(ierr, HYPRE_ERROR_CONV)) {
        HYPRE_GMRESGetMaxIter(sd->solver, &max_its);
        if (its <= max_its)
          cvg = CS_SLES_DIVERGED;
        else
          cvg = CS_SLES_MAX_ITERATION;
      }
    }
    break;

  case CS_SLES_HYPRE_FLEXGMRES:
    {
      /* Finalize setup and solve */
      HYPRE_FlexGMRESSetAbsoluteTol(sd->solver, precision*r_norm);
      HYPRE_FlexGMRESSetTol(sd->solver, 0);
      ierr = HYPRE_ParCSRFlexGMRESSolve(sd->solver, par_a, p_rhs, p_x);

      /* Get solution and information */
      HYPRE_ParCSRFlexGMRESGetFinalRelativeResidualNorm(sd->solver, &res);
      HYPRE_ParCSRFlexGMRESGetNumIterations(sd->solver, &its);

      /* Manage feedback */
      if (ierr != 0 && HYPRE_CheckError(ierr, HYPRE_ERROR_CONV)) {
        HYPRE_FlexGMRESGetMaxIter(sd->solver, &max_its);
        if (its <= max_its)
          cvg = CS_SLES_DIVERGED;
        else
          cvg = CS_SLES_MAX_ITERATION;
      }
    }
    break;

  case CS_SLES_HYPRE_LGMRES:
    {
      /* Finalize setup and solve */
      HYPRE_LGMRESSetAbsoluteTol(sd->solver, precision*r_norm);
      HYPRE_LGMRESSetTol(sd->solver, 0);
      ierr = HYPRE_ParCSRLGMRESSolve(sd->solver, par_a, p_rhs, p_x);

      /* Get solution and information */
      HYPRE_ParCSRLGMRESGetFinalRelativeResidualNorm(sd->solver, &res);
      HYPRE_ParCSRLGMRESGetNumIterations(sd->solver, &its);

      /* Manage feedback */
      if (ierr != 0 && HYPRE_CheckError(ierr, HYPRE_ERROR_CONV)) {
        HYPRE_LGMRESGetMaxIter(sd->solver, &max_its);
        if (its <= max_its)
          cvg = CS_SLES_DIVERGED;
        else
          cvg = CS_SLES_MAX_ITERATION;
      }
    }
    break;

  case CS_SLES_HYPRE_PCG:
    {
      /* Finalize setup and solve */
      HYPRE_PCGSetAbsoluteTol(sd->solver, precision*r_norm);
      HYPRE_PCGSetTol(sd->solver, 0);
      ierr = HYPRE_ParCSRPCGSolve(sd->solver, par_a, p_rhs, p_x);

      /* Get solution and information */
      HYPRE_ParCSRPCGGetFinalRelativeResidualNorm(sd->solver, &res);
      HYPRE_ParCSRPCGGetNumIterations(sd->solver, &its);

      /* Manage feedback */
      if (ierr != 0 && HYPRE_CheckError(ierr, HYPRE_ERROR_CONV)) {
        HYPRE_PCGGetMaxIter(sd->solver, &max_its);
        if (its <= max_its)
          cvg = CS_SLES_DIVERGED;
        else
          cvg = CS_SLES_MAX_ITERATION;
      }
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _("HYPRE: solver type (%s) not handled."),
              _cs_hypre_type_name(c->solver_type));
  }

  if (_t == NULL) {
    HYPRE_IJVectorGetValues(sd->coeffs->hx, n_rows, NULL, vx);
  }
  else {
    HYPRE_IJVectorGetValues(sd->coeffs->hx, n_rows, NULL, _t);
    for (HYPRE_BigInt ii = 0; ii < n_rows; ii++) {
      vx[ii] = _t[ii];
    }
    CS_FREE_HD(_t);
  }

  /* Set the convergence status. The feedback of HYPRE is quite minimal so that
     this only captures the convergence and avoids to get CS_SLES_ITERATING
     when an error has been raised. For some solvers, one can identify when the
     max. number of iterations has been reached. */
  if (ierr == 0)
    cvg = CS_SLES_CONVERGED;
  else {
    if (cvg == CS_SLES_ITERATING)
      cvg = CS_SLES_DIVERGED;
  }

  *residual = res;
  *n_iter = its;

  /* Update return values */
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
 * \brief Free HYPRE linear equation solver setup context.
 *
 * This function frees resolution-related data, such as
 * buffers and preconditioning but does not free the whole context,
 * as info used for logging (especially performance data) is maintained.
 *
 * \param[in, out]  context  pointer to iterative solver info and context
 *                           (actual type: cs_sles_hypre_t  *)
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_hypre_free(void  *context)
{
  cs_timer_t t0;
  t0 = cs_timer_time();

  cs_sles_hypre_t  *c  = context;

  if (c->setup_data != NULL) {

    cs_sles_hypre_setup_t *sd = c->setup_data;

    for (int i = 0; i < 2; i++) {
      HYPRE_Solver hs = (i == 0) ? sd->solver : sd->precond;
      cs_sles_hypre_type_t hs_type = (i == 0) ? c->solver_type : c->precond_type;

      if (hs == NULL)
        continue;

      switch(hs_type) {

      case CS_SLES_HYPRE_BOOMERAMG:
        HYPRE_BoomerAMGDestroy(hs);
        break;

      case CS_SLES_HYPRE_HYBRID:
        HYPRE_ParCSRHybridDestroy(hs);
        break;

      case CS_SLES_HYPRE_ILU:
        HYPRE_ILUDestroy(hs);
        break;

      case CS_SLES_HYPRE_BICGSTAB:
        HYPRE_ParCSRBiCGSTABDestroy(hs);
        break;

      case CS_SLES_HYPRE_GMRES:
        HYPRE_ParCSRGMRESDestroy(hs);
        break;

      case CS_SLES_HYPRE_FLEXGMRES:
        HYPRE_ParCSRFlexGMRESDestroy(hs);
        break;

      case CS_SLES_HYPRE_LGMRES:
        HYPRE_ParCSRLGMRESDestroy(hs);
        break;

      case CS_SLES_HYPRE_PCG:
        HYPRE_ParCSRPCGDestroy(hs);
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  _("HYPRE: solver type (%s) not handled."),
                  _cs_hypre_type_name(c->solver_type));
      }

      if (i == 0)
        sd->solver = NULL;
      else
        sd->precond = NULL;
    }

    BFT_FREE(c->setup_data);
  }

  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(c->t_setup), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the max. number of iterations associated to the given HYPRE
 *        contrext
 *
 * \param[in,out]  context       pointer to HYPRE linear solver info
 * \param[in]      n_max_iter    max. number of iterations
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_hypre_set_n_max_iter(cs_sles_hypre_t   *context,
                             int                n_max_iter)
{
  if (context == NULL)
    return;

  cs_sles_hypre_setup_t *sd = context->setup_data;

  if (sd == NULL)
    return; /* No need to continue. This will be done during the first call to
               the solve function */

  switch(context->solver_type) {

  case CS_SLES_HYPRE_BOOMERAMG:
    HYPRE_BoomerAMGSetMaxIter(sd->solver, n_max_iter);
    break;

  case CS_SLES_HYPRE_BICGSTAB:
    HYPRE_ParCSRBiCGSTABSetMaxIter(sd->solver, n_max_iter);
    break;

  case CS_SLES_HYPRE_GMRES:
    HYPRE_GMRESSetMaxIter(sd->solver, n_max_iter);
    break;

  case CS_SLES_HYPRE_FLEXGMRES:
    HYPRE_ParCSRFlexGMRESSetMaxIter(sd->solver, n_max_iter);
    break;

  case CS_SLES_HYPRE_LGMRES:
    HYPRE_ParCSRLGMRESSetMaxIter(sd->solver, n_max_iter);
    break;

  case CS_SLES_HYPRE_PCG:
    HYPRE_PCGSetMaxIter(sd->solver, n_max_iter);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Solver not handled yet", __func__);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define whether the solver should run on the host or
 *        accelerated device.
 *
 * If no device is available, this setting is ignored.
 *
 * \param[in,out]  context       pointer to HYPRE linear solver info
 * \param[in]      use_device    0 for host, 1 for device (GPU)
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_hypre_set_host_device(cs_sles_hypre_t   *context,
                              int                use_device)
{
  int device_id = cs_get_device_id();
  if (device_id < 0)
    use_device = 0;

  context->use_device = use_device;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Query whether the solver should run on the host or accelerated device.
 *
 * \param[in,out]  context       pointer to HYPRE linear solver info
 *
 * \return   0 for host, 1 for device (GPU)
 */
/*----------------------------------------------------------------------------*/

int
cs_sles_hypre_get_host_device(const cs_sles_hypre_t   *context)
{
  return context->use_device;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print information on hypre library.
 *
 * \param[in]  log_type  log type
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_hypre_library_info(cs_log_t  log_type)
{
 char hypre_config_options[256] = "";
 size_t l_max = 255;
 size_t l = l_max;

 snprintf(hypre_config_options, l, "%s %s (",
          HYPRE_RELEASE_NAME, HYPRE_RELEASE_VERSION);

 l -= strlen(hypre_config_options);

#if defined(HYPRE_USING_GPU)
 strncat(hypre_config_options, "GPU support, ", l);
 l -= strlen(hypre_config_options);
#endif

#if defined(HYPRE_SEQUENTIAL)
 strncat(hypre_config_options, "sequential, ", l);
 l -= strlen(hypre_config_options);
#endif

#if defined(HYPRE_USING_OPENMP)
 strncat(hypre_config_options, "OpenMP, ", l);
 l -= strlen(hypre_config_options);
#endif


#if defined(HYPRE_BIGINT)
 strncat(hypre_config_options, "large integers, ", l);
 l -= strlen(hypre_config_options);
#endif

 l = strlen(hypre_config_options) - 2;
 if (hypre_config_options[l] == ',') {
   hypre_config_options[l] = ')';
   hypre_config_options[l+1] = '\0';
 }
 else {
   hypre_config_options[l] = '\0';
 }

 cs_log_printf(log_type,
               "    %s\n", hypre_config_options);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
