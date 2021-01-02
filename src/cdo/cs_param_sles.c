/*============================================================================
 * Routines to handle the SLES settings
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2021 EDF S.A.

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
#include <stdlib.h>
#include <string.h>

#if defined(HAVE_PETSC)
#include <petscversion.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_error.h>
#include <bft_mem.h>

#include "cs_fp_exception.h"
#include "cs_log.h"
#include "cs_multigrid.h"
#include "cs_sles.h"

#if defined(HAVE_MUMPS)
#include "cs_sles_mumps.h"
#endif

#if defined(HAVE_PETSC)
#include "cs_sles_petsc.h"
#endif

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_param_sles.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Local private variables
 *============================================================================*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

#if defined(HAVE_PETSC)
/*----------------------------------------------------------------------------*/
/*!
 * \brief Predefined settings for GAMG as a preconditioner
 */
/*----------------------------------------------------------------------------*/

static inline void
_petsc_pcmg_hook(void)
{
#if PETSC_VERSION_GE(3,7,0)
  PetscOptionsSetValue(NULL, "-mg_levels_ksp_type", "richardson");
  PetscOptionsSetValue(NULL, "-mg_levels_pc_type", "sor");
  PetscOptionsSetValue(NULL, "-mg_levels_ksp_max_it", "1");
#else
  PetscOptionsSetValue("-mg_levels_ksp_type", "richardson");
  PetscOptionsSetValue("-mg_levels_pc_type", "sor");
  PetscOptionsSetValue("-mg_levels_ksp_max_it", "1");
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Predefined settings for GAMG as a preconditioner
 */
/*----------------------------------------------------------------------------*/

static inline void
_petsc_pcgamg_hook(void)
{
  _petsc_pcmg_hook();

#if PETSC_VERSION_GE(3,7,0)
  PetscOptionsSetValue(NULL, "-pc_gamg_threshold", "0.02");
  PetscOptionsSetValue(NULL, "-pc_gamg_reuse_interpolation", "TRUE");
  PetscOptionsSetValue(NULL, "-pc_gamg_square_graph", "4");
#else
  PetscOptionsSetValue("-pc_gamg_threshold", "0.02");
  PetscOptionsSetValue("-pc_gamg_reuse_interpolation", "TRUE");
  PetscOptionsSetValue("-pc_gamg_square_graph", "4");
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Predefined settings for BoomerAMG in HYPRE as a preconditioner
 */
/*----------------------------------------------------------------------------*/

static inline void
_petsc_pchypre_hook(void)
{
#if PETSC_VERSION_GE(3,7,0)
  PetscOptionsSetValue(NULL,"-pc_hypre_type","boomeramg");
  PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_coarsen_type","HMIS");
  PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_interp_type","ext+i-cc");
  PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_agg_nl","2");
  PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_P_max","4");
  PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_strong_threshold","0.5");
  PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_no_CF","");
#else
  PetscOptionsSetValue("-pc_hypre_type","boomeramg");
  PetscOptionsSetValue("-pc_hypre_boomeramg_coarsen_type","HMIS");
  PetscOptionsSetValue("-pc_hypre_boomeramg_interp_type","ext+i-cc");
  PetscOptionsSetValue("-pc_hypre_boomeramg_agg_nl","2");
  PetscOptionsSetValue("-pc_hypre_boomeramg_P_max","4");
  PetscOptionsSetValue("-pc_hypre_boomeramg_strong_threshold","0.5");
  PetscOptionsSetValue("-pc_hypre_boomeramg_no_CF","");
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set command line options for PC according to the kind of
 *        preconditionner
 *
 * \param[in]      slesp    set of parameters for the linear algebra
 * \param[in, out] pc       PETSc preconditioner structure
 */
/*----------------------------------------------------------------------------*/

static void
_petsc_set_pc_type(cs_param_sles_t   *slesp,
                   PC                 pc)
{
  if (slesp->solver == CS_PARAM_ITSOL_MUMPS ||
      slesp->solver == CS_PARAM_ITSOL_MUMPS_LDLT)
    return; /* Direct solver: Nothing to do at this stage */

  switch (slesp->precond) {

  case CS_PARAM_PRECOND_NONE:
    PCSetType(pc, PCNONE);
    break;

  case CS_PARAM_PRECOND_DIAG:
    PCSetType(pc, PCJACOBI);
    break;

  case CS_PARAM_PRECOND_BJACOB_ILU0:
  case CS_PARAM_PRECOND_BJACOB_SGS:
    PCSetType(pc, PCBJACOBI);
    break;

  case CS_PARAM_PRECOND_SSOR:
    PCSetType(pc, PCSOR);
    PCSORSetSymmetric(pc, SOR_SYMMETRIC_SWEEP);
    break;

  case CS_PARAM_PRECOND_ICC0:
#if defined(PETSC_HAVE_HYPRE)
    if (slesp->solver_class == CS_PARAM_SLES_CLASS_HYPRE) {
      PCSetType(pc, PCHYPRE);
      PCHYPRESetType(pc, "euclid");
    }
    else {
      PCSetType(pc, PCICC);
      PCFactorSetLevels(pc, 0);
    }
#else
    PCSetType(pc, PCICC);
    PCFactorSetLevels(pc, 0);
#endif
    break;

  case CS_PARAM_PRECOND_ILU0:
#if defined(PETSC_HAVE_HYPRE)
    if (slesp->solver_class == CS_PARAM_SLES_CLASS_HYPRE) {
      PCSetType(pc, PCHYPRE);
      PCHYPRESetType(pc, "euclid");
    }
    else {
      PCSetType(pc, PCBJACOBI);
    }
#else
    PCSetType(pc, PCBJACOBI);
#endif
    break;

  case CS_PARAM_PRECOND_AS:
    PCSetType(pc, PCASM);
    break;

  case CS_PARAM_PRECOND_AMG:
    {
      switch (slesp->amg_type) {

      case CS_PARAM_AMG_PETSC_GAMG:
        PCSetType(pc, PCGAMG);
        PCGAMGSetType(pc, PCGAMGAGG);
        PCGAMGSetNSmooths(pc, 1);
        break;

      case CS_PARAM_AMG_PETSC_PCMG:
        PCSetType(pc, PCMG);
        break;

      case CS_PARAM_AMG_HYPRE_BOOMER:
#if defined(PETSC_HAVE_HYPRE)
        PCSetType(pc, PCHYPRE);
        PCHYPRESetType(pc, "boomeramg");
#else
        cs_base_warn(__FILE__, __LINE__);
        cs_log_printf(CS_LOG_DEFAULT,
                      "%s: Eq. %s: Switch to MG since BoomerAMG is not"
                      " available.\n",
                      __func__, slesp->name);
#endif
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Eq. %s: Invalid AMG type for the PETSc library.",
                  __func__, slesp->name);
        break;

      } /* End of switch on the AMG type */

    } /* AMG as preconditioner */
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Eq. %s: Preconditioner not interfaced with PETSc.",
              __func__, slesp->name);
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set command line options for PC according to the kind of
 *        preconditionner
 *
 * \param[in]  slesp     set of parameters for the linear algebra solver
 */
/*----------------------------------------------------------------------------*/

static void
_petsc_set_pc_options_from_command_line(cs_param_sles_t   *slesp)
{
  switch (slesp->precond) {

#if defined(PETSC_HAVE_HYPRE)
  case CS_PARAM_PRECOND_ILU0:
  case CS_PARAM_PRECOND_ICC0:
    if (slesp->solver_class == CS_PARAM_SLES_CLASS_HYPRE) {
#if PETSC_VERSION_GE(3,7,0)
      PetscOptionsSetValue(NULL, "-pc_euclid_level", "-help");
#else
      PetscOptionsSetValue("-pc_euclid_level", "0");
#endif
    }
    break;
#endif  /* PETSc with HYPRE */

  case CS_PARAM_PRECOND_AMG:
    {
      switch (slesp->amg_type) {

      case CS_PARAM_AMG_PETSC_GAMG:
        _petsc_pcgamg_hook();
        break;

      case CS_PARAM_AMG_PETSC_PCMG:
        _petsc_pcmg_hook();
        break;

      case CS_PARAM_AMG_HYPRE_BOOMER:
#if defined(PETSC_HAVE_HYPRE)
        _petsc_pchypre_hook();
#else
        _petsc_pcmg_hook();
#endif
        break;

      default:
        break; /* Nothing else to do at this stage */

      } /* End of switch on the AMG type */

    } /* AMG as preconditioner */
    break;

  default:
    break; /* Nothing else to do at this stage */
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set PETSc solver
 *
 * \param[in]      slesp    pointer to SLES parameters
 * \param[in, out] ksp      pointer to PETSc KSP context
 */
/*----------------------------------------------------------------------------*/

static void
_petsc_set_krylov_solver(cs_param_sles_t    *slesp,
                         KSP                 ksp)
{
  /* 1) Set the type of normalization for the residual */
  switch (slesp->resnorm_type) {

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

  /* 2) Set the krylov solver */
  switch (slesp->solver) {

  case CS_PARAM_ITSOL_NONE:
    KSPSetType(ksp, KSPPREONLY);
    break;

  case CS_PARAM_ITSOL_BICG:      /* Improved Bi-CG stab */
    KSPSetType(ksp, KSPIBCGS);
    /* No choice otherwise PETSc yields an error */
    KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED);
    break;

  case CS_PARAM_ITSOL_BICGSTAB2: /* Preconditioned BiCGstab2 */
    KSPSetType(ksp, KSPBCGSL);
    break;

  case CS_PARAM_ITSOL_CG:        /* Preconditioned Conjugate Gradient */
    if (slesp->precond == CS_PARAM_PRECOND_AMG ||
        slesp->precond == CS_PARAM_PRECOND_AMG_BLOCK)
      KSPSetType(ksp, KSPFCG);
    else
      KSPSetType(ksp, KSPCG);
    break;

  case CS_PARAM_ITSOL_FCG:       /* Flexible Conjuguate Gradient */
    KSPSetType(ksp, KSPFCG);
    break;

  case CS_PARAM_ITSOL_FGMRES:     /* Preconditioned flexible GMRES */
    KSPSetType(ksp, KSPFGMRES);
    break;

  case CS_PARAM_ITSOL_GMRES:     /* Preconditioned GMRES */
    KSPSetType(ksp, KSPLGMRES);
    break;

  case CS_PARAM_ITSOL_MINRES:    /* Minimal residual */
    KSPSetType(ksp, KSPMINRES);
    break;

  case CS_PARAM_ITSOL_MUMPS:     /* Direct solver (factorization) */
  case CS_PARAM_ITSOL_MUMPS_LDLT:
#if defined(PETSC_HAVE_MUMPS)
    KSPSetType(ksp, KSPPREONLY);
#else
    bft_error(__FILE__, __LINE__, 0,
              " %s: MUMPS not interfaced with this installation of PETSc.",
              __func__);
#endif
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Iterative solver not interfaced with PETSc.", __func__);
  }

  /* 3) Additional settings arising from command lines */
  switch (slesp->solver) {

  case CS_PARAM_ITSOL_GMRES: /* Preconditioned GMRES */
#if PETSC_VERSION_GE(3,7,0)
    PetscOptionsSetValue(NULL, "-ksp_gmres_modifiedgramschmidt", "1");
#else
    PetscOptionsSetValue("-ksp_gmres_modifiedgramschmidt", "1");
#endif
    break;

  default:
    break; /* Nothing to do. Settings performed with another mechanism */

  }

  /* Apply modifications to the KSP structure given with command lines.
   * This setting stands for a first setting and may be overwritten with
   * parameters stored in the structure cs_param_sles_t
   *
   * Automatic monitoring
   *  PetscOptionsSetValue(NULL, "-ksp_monitor", "");
   *
   */

  KSPSetFromOptions(ksp);

  /* Apply settings from the cs_param_sles_t structure */
  switch (slesp->solver) {

  case CS_PARAM_ITSOL_GMRES: /* Preconditioned GMRES */
    {
      const int  n_max_restart = 40;

      KSPGMRESSetRestart(ksp, n_max_restart);
    }
    break;

#if defined(PETSC_HAVE_MUMPS)
  case CS_PARAM_ITSOL_MUMPS:
    {
      PC  pc;
      KSPGetPC(ksp, &pc);
      PCSetType(pc, PCLU);
      PCFactorSetMatSolverType(pc, MATSOLVERMUMPS);
    }
    break;

  case CS_PARAM_ITSOL_MUMPS_LDLT:
    {
      PC  pc;
      KSPGetPC(ksp, &pc);

      /* Retrieve the matrices related to this KSP */
      Mat a, pa;
      KSPGetOperators(ksp, &a, &pa);

      MatSetOption(a, MAT_SPD, PETSC_TRUE); /* set MUMPS id%SYM=1 */
      PCSetType(pc, PCCHOLESKY);

      PCFactorSetMatSolverType(pc, MATSOLVERMUMPS);
      PCFactorSetUpMatSolverType(pc); /* call MatGetFactor() to create F */
    }
    break;
#endif

  default:
    break; /* Nothing else to do */
  }

  /* Set KSP tolerances */
  PetscReal rtol, abstol, dtol;
  PetscInt  maxit;
  KSPGetTolerances(ksp, &rtol, &abstol, &dtol, &maxit);
  KSPSetTolerances(ksp,
                   slesp->eps,          /* relative convergence tolerance */
                   abstol,              /* absolute convergence tolerance */
                   dtol,                /* divergence tolerance */
                   slesp->n_max_iter);  /* max number of iterations */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set PETSc solver and preconditioner
 *
 * \param[in, out] context  pointer to optional (untyped) value or structure
 * \param[in, out] ksp      pointer to PETSc KSP context
 */
/*----------------------------------------------------------------------------*/

static void
_petsc_setup_hook(void   *context,
                  KSP     ksp)
{
  cs_param_sles_t  *slesp = (cs_param_sles_t  *)context;

  cs_fp_exception_disable_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */

  /* 1) Set the solver */
  _petsc_set_krylov_solver(slesp, ksp);

  /* Sanity checks */
  if (cs_glob_n_ranks > 1 && slesp->solver_class == CS_PARAM_SLES_CLASS_PETSC) {

    if (slesp->precond == CS_PARAM_PRECOND_ILU0 ||
        slesp->precond == CS_PARAM_PRECOND_ICC0) {
#if defined(PETSC_HAVE_HYPRE)
      slesp->solver_class = CS_PARAM_SLES_CLASS_HYPRE;
#else
      slesp->precond = CS_PARAM_PRECOND_BJACOB_ILU0;
      cs_base_warn(__FILE__, __LINE__);
      cs_log_printf(CS_LOG_DEFAULT,
                    " %s: System %s: Modify the requested preconditioner to"
                    " enable a parallel computation with PETSC.\n"
                    " Switch to a block jacobi preconditioner.\n"
                    " Please check your settings.", __func__, slesp->name);
#endif
    }
    else if (slesp->precond == CS_PARAM_PRECOND_SSOR) {

      slesp->precond = CS_PARAM_PRECOND_BJACOB_SGS;
      cs_base_warn(__FILE__, __LINE__);
      cs_log_printf(CS_LOG_DEFAULT,
                    " %s: System %s: Modify the requested preconditioner to"
                    " enable a parallel computation with PETSC.\n"
                    " Switch to a block jacobi preconditioner.\n"
                    " Please check your settings.", __func__, slesp->name);

    }

  } /* Advanced check for parallel run */

  /* 2) Set the preconditioner */
  PC  pc;
  KSPGetPC(ksp, &pc);
  _petsc_set_pc_type(slesp, pc);

  /* 3) Set PC options from command line */
  _petsc_set_pc_options_from_command_line(slesp);

  /* Apply modifications to the PC structure given with command lines.
   * This setting stands for a first setting and may be overwritten with
   * parameters stored in the structure cs_param_sles_t
   * To get the last word use cs_user_sles_petsc_hook()  */
  PCSetFromOptions(pc);

  /* 4) User function for additional settings */
  cs_user_sles_petsc_hook((void *)slesp, ksp);

  /* Dump the setup related to PETSc in a specific file */
  if (!slesp->setup_done) {

    KSPSetUp(ksp);
    cs_sles_petsc_log_setup(ksp);
    slesp->setup_done = true;

  }

  cs_fp_exception_restore_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Common settings for block preconditioning (when a system is split
 *         according to the Cartesian components: x,y,z)
 *
 * \param[in, out] slesp       pointer to the SLES parameter settings
 * \param[in, out] ksp         pointer to PETSc KSP structure (solver)
 * \param[in, out] pc          pointer to PETSc PC structure (preconditioner)
 * \param[in, out] xyz_subksp  list of KSP structures for each block
 */
/*----------------------------------------------------------------------------*/

static void
_petsc_common_block_hook(cs_param_sles_t    *slesp,
                         KSP                 ksp,
                         PC                  pc,
                         KSP               **xyz_subksp)
{
  /* Set KSP tolerances */
  PetscReal rtol, abstol, dtol;
  PetscInt  maxit;
  KSPGetTolerances(ksp, &rtol, &abstol, &dtol, &maxit);
  KSPSetTolerances(ksp,
                   slesp->eps,         /* relative convergence tolerance */
                   abstol,            /* absolute convergence tolerance */
                   dtol,              /* divergence tolerance */
                   slesp->n_max_iter); /* max number of iterations */

  switch (slesp->resnorm_type) {

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

  PCSetType(pc, PCFIELDSPLIT);
  PCFieldSplitSetType(pc, PC_COMPOSITE_ADDITIVE);

  /* Apply modifications to the KSP structure */
  PetscInt  id, n_split;

  PCFieldSplitSetBlockSize(pc, 3);
  id = 0;
  PCFieldSplitSetFields(pc, "x", 1, &id, &id);
  id = 1;
  PCFieldSplitSetFields(pc, "y", 1, &id, &id);
  id = 2;
  PCFieldSplitSetFields(pc, "z", 1, &id, &id);

  PCFieldSplitGetSubKSP(pc, &n_split, xyz_subksp);
  assert(n_split == 3);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer: setup hook for setting PETSc solver and
 *         preconditioner.
 *         Case of multiplicative AMG block preconditioner for a CG with GAMG
 *         as AMG type
 *
 * \param[in, out] context  pointer to optional (untyped) value or structure
 * \param[in, out] ksp      pointer to PETSc KSP context
 */
/*----------------------------------------------------------------------------*/

static void
_petsc_amg_block_gamg_hook(void     *context,
                           KSP       ksp)
{
  PC  pc;
  KSP  *xyz_subksp;
  cs_param_sles_t  *slesp = (cs_param_sles_t *)context;

  cs_fp_exception_disable_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */

  /* Set the solver */
  _petsc_set_krylov_solver(slesp, ksp);

  /* Common settings to block preconditionner */
  KSPGetPC(ksp, &pc);
  _petsc_common_block_hook(slesp,
                           ksp,
                           pc,
                           &xyz_subksp);

  /* Predefined settings when using AMG as a preconditioner */
  _petsc_pcgamg_hook();

  PC  _pc;
  for (PetscInt id = 0; id < 3; id++) {

    KSP  _ksp = xyz_subksp[id];
    KSPSetType(_ksp, KSPPREONLY);
    KSPGetPC(_ksp, &_pc);
    PCSetType(_pc, PCGAMG);

  }

  /* User function for additional settings */
  cs_user_sles_petsc_hook(context, ksp);

  PCSetFromOptions(pc);
  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);

  /* Dump the setup related to PETSc in a specific file */
  if (!slesp->setup_done) {
    cs_sles_petsc_log_setup(ksp);
    slesp->setup_done = true;
  }

  PetscFree(xyz_subksp);

  cs_fp_exception_restore_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */
}

#if defined(PETSC_HAVE_HYPRE)
/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer: setup hook for setting PETSc solver and
 *         preconditioner.
 *         Case of multiplicative AMG block preconditioner for a CG with boomer
 *         as AMG type
 *
 * \param[in, out] context  pointer to optional (untyped) value or structure
 * \param[in, out] ksp      pointer to PETSc KSP context
 */
/*----------------------------------------------------------------------------*/

static void
_petsc_amg_block_boomer_hook(void     *context,
                             KSP       ksp)
{
  cs_param_sles_t  *slesp = (cs_param_sles_t *)context;

  cs_fp_exception_disable_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */

  /* Set the solver */
  _petsc_set_krylov_solver(slesp, ksp);

  /* Common settings to block preconditionner */
  PC  pc;
  KSP  *xyz_subksp = NULL;

  KSPGetPC(ksp, &pc);
  _petsc_common_block_hook(slesp,
                           ksp,
                           pc,
                           &xyz_subksp);

  /* Predefined settings when using AMG as a preconditioner */
  _petsc_pchypre_hook();

  for (PetscInt id = 0; id < 3; id++) {

    KSP  _ksp = xyz_subksp[id];
    KSPSetType(_ksp, KSPPREONLY);

    PC  _pc;
    KSPGetPC(_ksp, &_pc);
    PCSetType(_pc, PCHYPRE);
    PCHYPRESetType(_pc, "boomeramg");

  }

  /* User function for additional settings */
  cs_user_sles_petsc_hook(context, ksp);

  PCSetFromOptions(pc);
  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);

  /* Dump the setup related to PETSc in a specific file */
  if (!slesp->setup_done) {
    cs_sles_petsc_log_setup(ksp);
    slesp->setup_done = true;
  }

  PetscFree(xyz_subksp);

  cs_fp_exception_restore_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */
}
#endif

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer: setup hook for setting PETSc solver and
 *         preconditioner.
 *         Case of block Jacobi preconditioner
 *
 * \param[in, out] context  pointer to optional (untyped) value or structure
 * \param[in, out] ksp      pointer to PETSc KSP context
 */
/*----------------------------------------------------------------------------*/

static void
_petsc_block_jacobi_hook(void     *context,
                         KSP       ksp)
{
  cs_param_sles_t  *slesp = (cs_param_sles_t *)context;

  cs_fp_exception_disable_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */

  /* Set the solver (tolerance and max_it too) */
  _petsc_set_krylov_solver(slesp, ksp);

  /* Common settings to block preconditionner */
  PC  pc;
  KSP  *xyz_subksp;

  KSPGetPC(ksp, &pc);
  _petsc_common_block_hook(slesp,
                           ksp,
                           pc,
                           &xyz_subksp);

  KSPSetUp(ksp);

  /* Predefined settings when using block-ILU as a preconditioner */
  PC  _pc;
  for (PetscInt id = 0; id < 3; id++) {

    KSP  _ksp = xyz_subksp[id];
    KSPSetType(_ksp, KSPPREONLY);
    KSPGetPC(_ksp, &_pc);
    if (slesp->solver_class == CS_PARAM_SLES_CLASS_HYPRE) {
#if defined(PETSC_HAVE_HYPRE)
      PCSetType(_pc, PCHYPRE);
      PCHYPRESetType(_pc, "euclid"); /* ILU(1) by default */
#else
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid option: HYPRE is not installed.", __func__);
#endif
    }
    else {
      PC  _subpc;
      KSP *_subksp;

      PCSetType(_pc, PCBJACOBI); /* Default for the block is an ILU(0) */
      KSPSetUp(_ksp);
      PCBJacobiGetSubKSP(_pc, NULL, NULL, &_subksp);
      KSPSetType(_subksp[0], KSPPREONLY);
      KSPGetPC(_subksp[0], &_subpc);

      if (slesp->precond == CS_PARAM_PRECOND_BJACOB_SGS)
        PCSetType(_subpc, PCEISENSTAT);
      else if (slesp->precond == CS_PARAM_PRECOND_BJACOB_ILU0) {
        PCFactorSetLevels(_pc, 0);
        PCFactorSetReuseOrdering(_pc, PETSC_TRUE);
        PCFactorSetMatOrderingType(_pc, MATORDERING1WD);
      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid preconditioner.",__func__);

    }

  }

  PCSetFromOptions(pc);
  PCSetUp(pc);

  /* User function for additional settings */
  cs_user_sles_petsc_hook(context, ksp);

  KSPSetFromOptions(ksp);

  /* Dump the setup related to PETSc in a specific file */
  if (!slesp->setup_done) {
    KSPSetUp(ksp);
    cs_sles_petsc_log_setup(ksp);
    slesp->setup_done = true;
  }

  PetscFree(xyz_subksp);

  cs_fp_exception_restore_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */
}
#endif /* defined(HAVE_PETSC) */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a \ref cs_param_sles_t structure and assign a default
 *         settings
 *
 * \param[in]  field_id      id related to to the variable field or -1
 * \param[in]  system_name   name of the system to solve or NULL
 *
 * \return a pointer to a cs_param_sles_t stucture
 */
/*----------------------------------------------------------------------------*/

cs_param_sles_t *
cs_param_sles_create(int          field_id,
                     const char  *system_name)
{
  cs_param_sles_t  *slesp = NULL;

  BFT_MALLOC(slesp, 1, cs_param_sles_t);

  slesp->verbosity = 0;                         /* SLES verbosity */
  slesp->field_id = field_id;                   /* associated field id */
  slesp->solver_class = CS_PARAM_SLES_CLASS_CS; /* solver family */
  slesp->precond = CS_PARAM_PRECOND_DIAG;       /* preconditioner */
  slesp->solver = CS_PARAM_ITSOL_GMRES;         /* iterative solver */
  slesp->amg_type = CS_PARAM_AMG_NONE;          /* no predefined AMG type */
  slesp->n_max_iter = 10000;                    /* max. number of iterations */
  slesp->eps = 1e-8;                            /* relative tolerance to stop
                                                   an iterative solver */
  slesp->resnorm_type = CS_PARAM_RESNORM_NONE;
  slesp->setup_done = false;

  slesp->name = NULL;
  if (system_name != NULL) {
    int  len = strlen(system_name) + 1;
    BFT_MALLOC(slesp->name, len, char);
    strncpy(slesp->name, system_name, len);
  }

  return slesp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a \ref cs_param_sles_t structure
 *
 * \param[in, out]  slesp    pointer to a \cs_param_sles_t structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_free(cs_param_sles_t   **p_slesp)
{
  if (p_slesp == NULL)
    return;

  cs_param_sles_t  *slesp = *p_slesp;

  if (slesp == NULL)
    return;

  BFT_FREE(slesp->name);
  BFT_FREE(slesp);
  slesp = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Log information related to the linear settings stored in the
 *         structure
 *
 * \param[in] slesp    pointer to a \ref cs_param_sles_log
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_log(cs_param_sles_t   *slesp)
{
  if (slesp == NULL)
    return;

  cs_log_printf(CS_LOG_SETUP, "\n### %s | Linear algebra settings\n",
                slesp->name);
  cs_log_printf(CS_LOG_SETUP, "  * %s | SLES Family:", slesp->name);
  if (slesp->solver_class == CS_PARAM_SLES_CLASS_CS)
    cs_log_printf(CS_LOG_SETUP, "             Code_Saturne\n");
  else if (slesp->solver_class == CS_PARAM_SLES_CLASS_MUMPS)
    cs_log_printf(CS_LOG_SETUP, "             MUMPS\n");
  else if (slesp->solver_class == CS_PARAM_SLES_CLASS_HYPRE)
    cs_log_printf(CS_LOG_SETUP, "             HYPRE\n");
  else if (slesp->solver_class == CS_PARAM_SLES_CLASS_PETSC)
    cs_log_printf(CS_LOG_SETUP, "             PETSc\n");

  cs_log_printf(CS_LOG_SETUP, "  * %s | SLES Verbosity:          %d\n",
                slesp->name, slesp->verbosity);
  cs_log_printf(CS_LOG_SETUP, "  * %s | SLES Field id:           %d\n",
                slesp->name, slesp->field_id);
  cs_log_printf(CS_LOG_SETUP, "  * %s | SLES Solver.MaxIter:     %d\n",
                slesp->name, slesp->n_max_iter);

  cs_log_printf(CS_LOG_SETUP, "  * %s | SLES Solver.Name:        %s\n",
                slesp->name, cs_param_get_solver_name(slesp->solver));
  if (slesp->solver == CS_PARAM_ITSOL_AMG)
    cs_log_printf(CS_LOG_SETUP, "  * %s | SLES AMG.Type:           %s\n",
                  slesp->name, cs_param_get_amg_type_name(slesp->amg_type));
  cs_log_printf(CS_LOG_SETUP, "  * %s | SLES Solver.Precond:     %s\n",
                slesp->name, cs_param_get_precond_name(slesp->precond));
  if (slesp->precond == CS_PARAM_PRECOND_AMG)
    cs_log_printf(CS_LOG_SETUP, "  * %s | SLES AMG.Type:           %s\n",
                  slesp->name, cs_param_get_amg_type_name(slesp->amg_type));

  cs_log_printf(CS_LOG_SETUP, "  * %s | SLES Solver.Eps:        % -10.6e\n",
                slesp->name, slesp->eps);

  cs_log_printf(CS_LOG_SETUP, "  * %s | SLES Normalization:      ",
                slesp->name);
  switch (slesp->resnorm_type) {
  case CS_PARAM_RESNORM_NORM2_RHS:
    cs_log_printf(CS_LOG_SETUP, "Euclidean norm of the RHS\n");
    break;
  case CS_PARAM_RESNORM_WEIGHTED_RHS:
    cs_log_printf(CS_LOG_SETUP, "Weighted Euclidean norm of the RHS\n");
    break;
  case CS_PARAM_RESNORM_FILTERED_RHS:
    cs_log_printf(CS_LOG_SETUP, "Filtered Euclidean norm of the RHS\n");
    break;
  case CS_PARAM_RESNORM_NONE:
  default:
    cs_log_printf(CS_LOG_SETUP, "None\n");
    break;
  }
  cs_log_printf(CS_LOG_SETUP, "\n");
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Copy a cs_param_sles_t structure from src to dst
 *
 * \param[in]       src    reference cs_param_sles_t structure to copy
 * \param[in, out]  dst    copy of the reference at exit
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_copy_from(cs_param_sles_t   *src,
                        cs_param_sles_t   *dst)
{
  if (dst == NULL)
    return;

  /* Remark: name is managed at the creation of the structure */

  dst->setup_done = src->setup_done;
  dst->verbosity = src->verbosity;
  dst->field_id = src->field_id;

  dst->solver_class = src->solver_class;
  dst->precond = src->precond;
  dst->solver = src->solver;
  dst->amg_type = src->amg_type;

  dst->resnorm_type = src->resnorm_type;
  dst->n_max_iter = src->n_max_iter;
  dst->eps = src->eps;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set parameters for initializing SLES structures used for the
 *        resolution of the linear system.
 *        Case of saturne's own solvers.
 *
 * \param[in]       use_field_id  if false use system name
 * \param[in, out]  slesp         pointer to a \ref cs_param_sles_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_param_set_saturne_sles(bool                 use_field_id,
                                   cs_param_sles_t     *slesp)
{
  assert(slesp != NULL);  /* Sanity checks */

  const  char  *sles_name = use_field_id ? NULL : slesp->name;
  assert(slesp->field_id > -1 || sles_name != NULL);

  /* 1- Define the preconditioner */
  /*    ========================= */

  int  poly_degree;
  cs_sles_pc_t *pc = NULL;

  switch (slesp->precond) {

  case CS_PARAM_PRECOND_DIAG:
    poly_degree = 0;
    break;

  case CS_PARAM_PRECOND_POLY1:
    poly_degree = 1;
    break;

  case CS_PARAM_PRECOND_POLY2:
    poly_degree = 2;
    break;

  case CS_PARAM_PRECOND_AMG:
    poly_degree = -1;
    switch (slesp->amg_type) {

    case CS_PARAM_AMG_HOUSE_V:
      pc = cs_multigrid_pc_create(CS_MULTIGRID_V_CYCLE);
      break;
    case CS_PARAM_AMG_HOUSE_K:
      if (slesp->solver == CS_PARAM_ITSOL_CG)
        slesp->solver = CS_PARAM_ITSOL_FCG;
      pc = cs_multigrid_pc_create(CS_MULTIGRID_K_CYCLE);
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " %s: System: %s; Invalid AMG type with Code_Saturne solvers.",
                __func__, slesp->name);
      break;

    } /* Switch on the type of AMG */
    break; /* AMG as preconditioner */

  case CS_PARAM_PRECOND_GKB_CG:
  case CS_PARAM_PRECOND_GKB_GMRES:
    poly_degree = -1;
    break;

  case CS_PARAM_PRECOND_NONE:
  default:
    poly_degree = -1;       /* None or other */

  } /* Switch on the preconditioner type */

  /* 2- Define the iterative solver */
  /*    =========================== */

  cs_sles_it_t  *it = NULL;
  cs_multigrid_t  *mg = NULL;

  switch (slesp->solver) {

  case CS_PARAM_ITSOL_AMG:
    {
      switch (slesp->amg_type) {

      case CS_PARAM_AMG_HOUSE_V:
        mg = cs_multigrid_define(slesp->field_id,
                                 sles_name,
                                 CS_MULTIGRID_V_CYCLE);

        /* Advanced setup (default is specified inside the brackets)
         * for AMG as solver */
        cs_multigrid_set_solver_options
          (mg,
           CS_SLES_JACOBI,   /* descent smoother type (CS_SLES_PCG) */
           CS_SLES_JACOBI,   /* ascent smoother type (CS_SLES_PCG) */
           CS_SLES_PCG,      /* coarse solver type (CS_SLES_PCG) */
           slesp->n_max_iter, /* n max cycles (100) */
           5,                /* n max iter for descent (10) */
           5,                /* n max iter for ascent (10) */
           1000,             /* n max iter coarse solver (10000) */
           0,                /* polynomial precond. degree descent (0) */
           0,                /* polynomial precond. degree ascent (0) */
           -1,               /* polynomial precond. degree coarse (0) */
           1.0,    /* precision multiplier descent (< 0 forces max iters) */
           1.0,    /* precision multiplier ascent (< 0 forces max iters) */
           1);     /* requested precision multiplier coarse (default 1) */
        break;

      case CS_PARAM_AMG_HOUSE_K:
        mg = cs_multigrid_define(slesp->field_id,
                                 sles_name,
                                 CS_MULTIGRID_K_CYCLE);

        cs_multigrid_set_solver_options
          (mg,
           CS_SLES_P_SYM_GAUSS_SEIDEL, /* descent smoother */
           CS_SLES_P_SYM_GAUSS_SEIDEL, /* ascent smoother */
           CS_SLES_PCG,                /* coarse smoother */
           slesp->n_max_iter,           /* n_max_cycles */
           1,                          /* n_max_iter_descent, */
           1,                          /* n_max_iter_ascent */
           100,                        /* n_max_iter_coarse */
           0,                          /* poly_degree_descent */
           0,                          /* poly_degree_ascent */
           0,                          /* poly_degree_coarse */
           -1.0,                       /* precision_mult_descent */
           -1.0,                       /* precision_mult_ascent */
           1);                         /* precision_mult_coarse */
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  " %s; System: %s -- Invalid AMG type with Code_Saturne"
                  " solvers.", __func__, slesp->name);
        break;

      } /* End of switch on the AMG type */

    }
    break; /* AMG as solver */

  case CS_PARAM_ITSOL_BICG:
    it = cs_sles_it_define(slesp->field_id,
                           sles_name,
                           CS_SLES_BICGSTAB,
                           poly_degree,
                           slesp->n_max_iter);
    break;

  case CS_PARAM_ITSOL_BICGSTAB2:
    it = cs_sles_it_define(slesp->field_id,
                           sles_name,
                           CS_SLES_BICGSTAB2,
                           poly_degree,
                           slesp->n_max_iter);
    break;

  case CS_PARAM_ITSOL_CG:
    it = cs_sles_it_define(slesp->field_id,
                           sles_name,
                           CS_SLES_PCG,
                           poly_degree,
                           slesp->n_max_iter);
    break;

  case CS_PARAM_ITSOL_CR3:
    it = cs_sles_it_define(slesp->field_id,
                           sles_name,
                           CS_SLES_PCR3,
                           poly_degree,
                           slesp->n_max_iter);
    break;

  case CS_PARAM_ITSOL_FCG:
    it = cs_sles_it_define(slesp->field_id,
                           sles_name,
                           CS_SLES_IPCG,
                           poly_degree,
                           slesp->n_max_iter);
    break;

  case CS_PARAM_ITSOL_GAUSS_SEIDEL:
    it = cs_sles_it_define(slesp->field_id,
                           sles_name,
                           CS_SLES_P_GAUSS_SEIDEL,
                           -1, /* Not useful to apply a preconditioner */
                           slesp->n_max_iter);
    break;

  case CS_PARAM_ITSOL_GKB_CG:
    it = cs_sles_it_define(slesp->field_id,
                           sles_name,
                           CS_SLES_IPCG, /* Flexible CG */
                           poly_degree,
                           slesp->n_max_iter);
    break;

  case CS_PARAM_ITSOL_GKB_GMRES:
    it = cs_sles_it_define(slesp->field_id,
                           sles_name,
                           CS_SLES_GMRES, /* Should a flexible GMRES */
                           poly_degree,
                           slesp->n_max_iter);
    break;

  case CS_PARAM_ITSOL_GMRES:
    it = cs_sles_it_define(slesp->field_id,
                           sles_name,
                           CS_SLES_GMRES,
                           poly_degree,
                           slesp->n_max_iter);
    break;

  case CS_PARAM_ITSOL_JACOBI:
    it = cs_sles_it_define(slesp->field_id,
                           sles_name,
                           CS_SLES_JACOBI,
                           -1, /* Not useful to apply a preconditioner */
                           slesp->n_max_iter);
    break;

  case CS_PARAM_ITSOL_SYM_GAUSS_SEIDEL:
    it = cs_sles_it_define(slesp->field_id,
                           sles_name,
                           CS_SLES_P_SYM_GAUSS_SEIDEL,
                           -1, /* Not useful to apply a preconditioner */
                           slesp->n_max_iter);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid iterative solver for solving equation %s.\n"
              " Please modify your settings.", __func__, slesp->name);
    break;

  } /* end of switch */

  /* Update the preconditioner settings if needed */
  if (slesp->precond == CS_PARAM_PRECOND_AMG) {

    assert(pc != NULL && it != NULL);

    mg = cs_sles_pc_get_context(pc);
    cs_sles_it_transfer_pc(it, &pc);

    /* Change the default settings for CDO/HHO when used as preconditioner */
    cs_multigrid_set_solver_options
      (mg,
       CS_SLES_PCG,       /* descent smoother */
       CS_SLES_PCG,       /* ascent smoother */
       CS_SLES_PCG,       /* coarse solver */
       slesp->n_max_iter, /* n_max_cycles */
       4,                 /* n_max_iter_descent, */
       4,                 /* n_max_iter_ascent */
       200,               /* n_max_iter_coarse */
       0,                 /* poly_degree_descent */
       0,                 /* poly_degree_ascent */
       0,                 /* poly_degree_coarse */
       -1.0,              /* precision_mult_descent */
       -1.0,              /* precision_mult_ascent */
       1.0);              /* precision_mult_coarse */

    /* If this is a K-cycle multigrid. Change the default aggregation
       algorithm */
    if (slesp->amg_type == CS_PARAM_AMG_HOUSE_K)
      cs_multigrid_set_coarsening_options(mg,
                                          8,   /* aggregation_limit*/
                                          CS_GRID_COARSENING_SPD_MX,
                                          10,  /* n_max_levels */
                                          50,  /* min_g_cells */
                                          0.,  /* P0P1 relaxation */
                                          0);  /* postprocess */

  } /* AMG as preconditioner */

  /* Define the level of verbosity for SLES structure */
  if (slesp->verbosity > 3) {

    cs_sles_t  *sles = cs_sles_find_or_add(slesp->field_id, sles_name);
    cs_sles_it_t  *sles_it = (cs_sles_it_t *)cs_sles_get_context(sles);

    /* true = use_iteration instead of wall clock time */
    cs_sles_it_set_plot_options(sles_it, slesp->name, true);

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set parameters for initializing SLES structures used for the
 *        resolution of the linear system.
 *        Case of MUMPS's own solvers.
 *
 * \param[in]       use_field_id  if false use system name
 * \param[in, out]  slesp         pointer to a \ref cs_param_sles_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_param_set_mumps_sles(bool                 use_field_id,
                                 cs_param_sles_t     *slesp)
{
  assert(slesp != NULL);  /* Sanity checks */

  const  char  *sles_name = use_field_id ? NULL : slesp->name;
  assert(slesp->field_id > -1 || sles_name != NULL);

  int  sym = 0; /* One assumes a non-symmetric as default */
  if (slesp->solver == CS_PARAM_ITSOL_MUMPS_LDLT)
    sym = 1;

#if defined(HAVE_MUMPS)
  cs_sles_mumps_define(slesp->field_id,
                       sles_name,
                       sym,
                       slesp->verbosity,
                       cs_user_sles_mumps_hook,
                       (void *)slesp);
#else
  bft_error(__FILE__, __LINE__, 0,
            "%s: System: %s\n"
            " MUMPS is not supported directly.\n"
            " Please check your settings or your code_saturne installation.",
            __func__, slesp->name);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set parameters for initializing SLES structures used for the
 *        resolution of the linear system.
 *        Case of PETSc and Hypre families of solvers.
 *
 * \param[in]       use_field_id  if false use system name
 * \param[in, out]  slesp         pointer to a \ref cs_param_sles_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_param_set_petsc_hypre_sles(bool                 use_field_id,
                                       cs_param_sles_t     *slesp)
{
  assert(slesp != NULL);  /* Sanity checks */

  const  char  *sles_name = use_field_id ? NULL : slesp->name;
  assert(slesp->field_id > -1 || sles_name != NULL);

#if defined(HAVE_PETSC)
  cs_sles_petsc_init();

  if (slesp->precond == CS_PARAM_PRECOND_AMG_BLOCK) {

    if (slesp->amg_type == CS_PARAM_AMG_PETSC_GAMG) {
      cs_sles_petsc_define(slesp->field_id,
                           sles_name,
                           MATMPIAIJ,
                           _petsc_amg_block_gamg_hook,
                           (void *)slesp);
    }
    else if (slesp->amg_type == CS_PARAM_AMG_HYPRE_BOOMER) {
#if defined(PETSC_HAVE_HYPRE)
      cs_sles_petsc_define(slesp->field_id,
                           sles_name,
                           MATMPIAIJ,
                           _petsc_amg_block_boomer_hook,
                           (void *)slesp);
#else
      bft_error(__FILE__, __LINE__, 0,
                " %s: System: %s.\n"
                " Boomer is not available. Switch to another AMG solver.",
                __func__, slesp->name);
#endif
    }
    else
      bft_error(__FILE__, __LINE__, 0,
                " %s: System: %s\n"
                " No AMG solver available for a block-AMG.", __func__,
                slesp->name);

  }
  else if (slesp->precond == CS_PARAM_PRECOND_BJACOB_ILU0 ||
           slesp->precond == CS_PARAM_PRECOND_BJACOB_SGS)
    cs_sles_petsc_define(slesp->field_id,
                         sles_name,
                         MATMPIAIJ,
                         _petsc_block_jacobi_hook,
                         (void *)slesp);
  else
    cs_sles_petsc_define(slesp->field_id,
                         sles_name,
                         MATMPIAIJ,
                         _petsc_setup_hook,
                         (void *)slesp);

#else
    bft_error(__FILE__, __LINE__, 0,
              _(" %s: PETSC algorithms used to solve %s are not linked.\n"
                " Please install Code_Saturne with PETSc."),
              __func__, slesp->name);
#endif /* HAVE_PETSC */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define cs_sles_t structure in accordance with the settings of a
 *        cs_param_sles_t structure (SLES = Sparse Linear Equation Solver)
 *
 * \param[in]       use_field_id  if false use system name to define a SLES
 * \param[in, out]  slesp         pointer to a cs_param_sles_t structure
 *
 * \return an error code (-1 if a problem is encountered, 0 otherwise)
 */
/*----------------------------------------------------------------------------*/

int
cs_param_sles_set(bool                 use_field_id,
                  cs_param_sles_t     *slesp)
{
  if (slesp == NULL)
    return 0;

  switch (slesp->solver_class) {

  case CS_PARAM_SLES_CLASS_CS: /* Code_Saturne's own solvers */
    /* true = use field_id instead of slesp->name to set the sles */
    cs_equation_param_set_saturne_sles(use_field_id, slesp);
    break;

  case CS_PARAM_SLES_CLASS_MUMPS: /* MUMPS sparse direct solvers */
    /* true = use field_id instead of slesp->name to set the sles */
    cs_equation_param_set_mumps_sles(use_field_id, slesp);
    break;

  case CS_PARAM_SLES_CLASS_PETSC: /* PETSc solvers */
  case CS_PARAM_SLES_CLASS_HYPRE: /* HYPRE solvers through PETSc */
    /* true = use field_id instead of slesp->name to set the sles */
    cs_equation_param_set_petsc_hypre_sles(use_field_id, slesp);
    break;

  default:
    return -1;

  } /* End of switch on class of solvers */

  /* Define the level of verbosity for the SLES structure */
  if (slesp->verbosity > 1) {

    /* All the previous SLES are defined thanks to the field_id */
    cs_sles_t  *sles = NULL;
    if (use_field_id)
      sles = cs_sles_find_or_add(slesp->field_id, NULL);
    else
      sles = cs_sles_find_or_add(slesp->field_id, slesp->name);

    /* Set verbosity */
    cs_sles_set_verbosity(sles, slesp->verbosity);

  }

  return 0;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
