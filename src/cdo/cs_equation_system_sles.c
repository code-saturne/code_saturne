/*============================================================================
 * Functions to handle the solving step of systems of equations hinging on the
 * cs_equation_system_t structure
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

#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_equation_param.h"
#include "cs_fp_exception.h"
#include "cs_param_sles.h"

#if defined(HAVE_MUMPS)
#include "cs_sles_mumps.h"
#endif

#if defined(HAVE_PETSC)
#include "cs_sles_petsc.h"
#endif

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_equation_system_sles.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local macro definitions
 *============================================================================*/

#define CS_EQUATION_SYSTEM_SLES_DBG  0

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Local variables
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private variables
 *============================================================================*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

#if defined(HAVE_PETSC)
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
  cs_equation_system_param_t  *sysp = context;
  cs_iter_algo_param_t  sp = sysp->linear_solver;

  cs_fp_exception_disable_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */

  PC  pc;

  KSPSetType(ksp, KSPPREONLY);
  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCLU);
  PCFactorSetMatSolverType(pc, MATSOLVERMUMPS);

  KSPSetTolerances(ksp,
                   sp.rtol,             /* relative convergence tolerance */
                   sp.atol,             /* absolute convergence tolerance */
                   sp.dtol,             /* divergence tolerance */
                   sp.n_max_algo_iter); /* max number of iterations */

  /* Dump the setup related to PETSc in a specific file */

  if (!sysp->sles_setup_done) {
    cs_sles_petsc_log_setup(ksp);
    sysp->sles_setup_done = true;
  }

  cs_fp_exception_restore_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */
}
#endif  /* PETSC_HAVE_MUMPS */
#endif  /* HAVE_PETSC */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the SLES associated to a system of equation
 *
 * \param[in]  n_eqs     number of equations in the system to solve
 * \param[in]  sysp      set of paremeters for the system of equations
 * \param[in]  blocks    array of the core members for an equation
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_system_sles_init(int                            n_eqs,
                             cs_equation_system_param_t    *sysp,
                             cs_equation_core_t           **blocks)
{
  CS_UNUSED(n_eqs);
  assert(sysp != NULL);

  switch (sysp->sles_strategy) {

  case CS_EQUATION_SYSTEM_SLES_MUMPS:
    {
#if defined(HAVE_MUMPS)
      cs_equation_core_t  *block00 = blocks[0];
      cs_equation_param_t  *eqp00 = block00->param;
      cs_param_sles_t  *slesp00 = eqp00->sles_param;

      if (slesp00->solver != CS_PARAM_ITSOL_MUMPS &&
          slesp00->solver != CS_PARAM_ITSOL_MUMPS_LDLT &&
          slesp00->solver != CS_PARAM_ITSOL_MUMPS_FLOAT &&
          slesp00->solver != CS_PARAM_ITSOL_MUMPS_FLOAT_LDLT)
        slesp00->solver = CS_PARAM_ITSOL_MUMPS;

      cs_sles_mumps_define(-1,
                           sysp->name,
                           slesp00,
                           cs_user_sles_mumps_hook,
                           NULL);
#else
#if defined(HAVE_PETSC)
#if defined(PETSC_HAVE_MUMPS)
      cs_sles_petsc_init();
      cs_sles_petsc_define(-1,
                           sysp->name,
                           MATMPIAIJ,
                           _mumps_hook,
                           (void *)sysp);
#else
      bft_error(__FILE__, __LINE__, 0,
                "%s: Invalid strategy for solving the linear system %s\n"
                " PETSc with MUMPS is required with this option.\n",
                __func__, sysp->name);
#endif  /* PETSC_HAVE_MUMPS */
      bft_error(__FILE__, __LINE__, 0,
                "%s: Invalid strategy for solving the linear system %s\n"
                " Neither PETSc nor MUMPS is available.\n",
                __func__, sysp->name);
#endif  /* HAVE_PETSC */
#endif  /* HAVE_MUMPS */
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid strategy for solving the linear system \"%s\"\n",
              __func__, sysp->name);

  } /* switch on the SLES strategy */

  sysp->sles_setup_done = true;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
