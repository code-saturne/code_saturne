/*============================================================================
 * Routines to handle the SLES settings for a saddle-point problem
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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
#include <float.h>

#if defined(HAVE_PETSC)
#include <petsc.h>
#include <petscconf.h> /* Useful to know if HYPRE is accessible through PETSc */
#include <petscversion.h>
#endif

/* Avoid warnings due to previous values */
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_URL
#undef PACKAGE_VERSION

#if defined(HAVE_HYPRE)
#include <HYPRE_krylov.h>
#include <HYPRE_parcsr_ls.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_error.h>
#include <bft_mem.h>

#include "cs_equation.h"
#include "cs_fp_exception.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_multigrid.h"
#include "cs_param_sles_setup.h"
#include "cs_saddle_solver.h"
#include "cs_sles.h"

#if defined(HAVE_MUMPS)
#include "cs_sles_mumps.h"
#endif

#if defined(HAVE_PETSC)
#include "cs_sles_petsc.h"
#endif

#if defined(HAVE_HYPRE)
#include "cs_sles_hypre.h"
#endif

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_saddle_solver_setup.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

#if defined(HAVE_PETSC)
/*----------------------------------------------------------------------------*/
/*!
 * \brief Generate two IndexSet structures for the definition of a FieldSplit
 *        preconditioner
 *
 * \param[in]      solver  pointer to a saddle-point solver
 * \param[in, out] is1     IndexSet for the (1,1)-block DoFs
 * \param[in, out] is2     IndexSet for the (2,2)-block DoFs
 */
/*----------------------------------------------------------------------------*/

static void
_build_is_for_fieldsplit(cs_saddle_solver_t  *solver,
                         IS                  *is1,
                         IS                  *is2)
{
  cs_cdo_system_helper_t  *sh = solver->system_helper;

  const cs_lnum_t  n1_dofs = solver->n1_scatter_dofs;
  const cs_lnum_t  n2_dofs = solver->n2_scatter_dofs;

  /* The block id equal to 0 should describe the full system */

  const cs_range_set_t  *rset = cs_cdo_system_get_range_set(sh, 0);

  PetscInt  *indices = nullptr;
  PetscMalloc1(CS_MAX(n1_dofs, n2_dofs), &indices);

  /* IndexSet for the (1,1)-block of DoFs */

  if (rset->n_elts[0] == rset->n_elts[1]) { /* Keep all DoFs */

    for (PetscInt i = 0; i < n1_dofs; i++)
      indices[i] = rset->g_id[i];

    ISCreateGeneral(PETSC_COMM_SELF, n1_dofs, indices, PETSC_COPY_VALUES, is1);

  }
  else { /* Filter DoFs */

    PetscInt  n1_keep = 0;
    for (PetscInt i = 0; i < n1_dofs; i++) {

      cs_gnum_t  g_id = rset->g_id[i];
      if (g_id >= rset->l_range[0] && g_id < rset->l_range[1])
        indices[n1_keep++] = g_id;

    }

    ISCreateGeneral(PETSC_COMM_WORLD, n1_keep, indices, PETSC_COPY_VALUES, is1);

  }

#if defined(DEBUG) && !defined(NDEBUG) && CS_SADDLE_SOLVER_SETUP_DBG > 1
  /* Print the index set to stdout */
  ISView(*is1, PETSC_VIEWER_STDOUT_SELF);
#endif

  /* Re-used the buffer for indices to create the second IndexSet The second
   * set of DoFs is assumed to be located in cells. This is true for all
   * discretization so far: CDO-Fb schemes for NavSto, CDO-Cb schemes for
   * scalar and MAC schemes. So, the treatment should be the same in sequential
   * and parallel computation
   */

  const cs_gnum_t  *g2_ids = rset->g_id + n1_dofs;
  for (PetscInt i = 0; i < n2_dofs; i++)
    indices[i] = g2_ids[i];

  ISCreateGeneral(PETSC_COMM_SELF, n2_dofs, indices, PETSC_COPY_VALUES, is2);

#if defined(DEBUG) && !defined(NDEBUG) && CS_SADDLE_SOLVER_SETUP_DBG > 1
  /* Print the index set to stdout */
  ISView(*is2, PETSC_VIEWER_STDOUT_SELF);
#endif

  PetscFree(indices);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the corresponding PC associated to the given set of parameters
 *
 * \param[in] slesp  set of parameters for the linear algebra
 *
 * \return a type of PETSc preconditioner
 */
/*----------------------------------------------------------------------------*/

/* static PCType */
/* _petsc_get_pc_type(const cs_param_sles_t  *slesp) */
/* { */
/*   PCType  pc_type = PCNONE; */

/*   switch (slesp->precond) { */

/*   case CS_PARAM_PRECOND_NONE: */
/*     return PCNONE; */

/*   case CS_PARAM_PRECOND_DIAG: */
/*     return PCJACOBI; */

/*   case CS_PARAM_PRECOND_BJACOB_ILU0: */
/*   case CS_PARAM_PRECOND_BJACOB_SGS: */
/*     return PCBJACOBI; */

/*   case CS_PARAM_PRECOND_SSOR: */
/*     return PCSOR; */

/*   case CS_PARAM_PRECOND_ICC0: */
/*     return PCICC; */

/*   case CS_PARAM_PRECOND_ILU0: */
/*     return PCILU; */

/*   case CS_PARAM_PRECOND_LU: */
/*     return PCLU; */

/*   case CS_PARAM_PRECOND_AMG: */
/*     { */
/*       switch (slesp->amg_type) { */

/*       case CS_PARAM_AMG_PETSC_GAMG_V: */
/*       case CS_PARAM_AMG_PETSC_GAMG_W: */
/*         return PCGAMG; */
/*         break; */

/*       case CS_PARAM_AMG_PETSC_PCMG: */
/*         return PCMG; */
/*         break; */

/*       case CS_PARAM_AMG_HYPRE_BOOMER_V: */
/*       case CS_PARAM_AMG_HYPRE_BOOMER_W: */
/*         if (cs_param_sles_hypre_from_petsc()) */
/*           return PCHYPRE; */

/*         else { */
/*           cs_base_warn(__FILE__, __LINE__); */
/*           cs_log_printf(CS_LOG_DEFAULT, */
/*                         "%s: Switch to MG since BoomerAMG is not available.\n", */
/*                         __func__); */
/*           return PCMG; */
/*         } */
/*         break; */

/*       default: */
/*         bft_error(__FILE__, __LINE__, 0, */
/*                   " %s: Invalid AMG type for the PETSc library.", __func__); */
/*         break; */

/*       } /\* End of switch on the AMG type *\/ */

/*     } /\* AMG as preconditioner *\/ */
/*     break; */

/*   default: */
/*     bft_error(__FILE__, __LINE__, 0, */
/*               " %s: Preconditioner not interfaced with PETSc.", __func__); */
/*     break; */
/*   } */

/*   return pc_type; */
/* } */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function pointer: setup hook for setting PETSc solver and
 *        preconditioner.
 *
 * \param[in]      slesp  pointer to a set of SLES settings
 * \param[in, out] u_ksp  pointer to PETSc KSP context
 */
/*----------------------------------------------------------------------------*/

/* static void */
/* _set_vector_block11_ksp(cs_param_sles_t  *slesp, */
/*                         KSP               u_ksp) */
/* { */
/*   PC  u_pc; */
/*   KSPGetPC(u_ksp, &u_pc); */
/*   PCType  pc_type = _petsc_get_pc_type(slesp); */

/*   /\* Set the type of solver *\/ */

/*   switch (slesp->solver) { */

/*   case CS_PARAM_ITSOL_NONE: */
/*     KSPSetType(u_ksp, KSPPREONLY); */
/*     break; */
/*   case CS_PARAM_ITSOL_FCG: */
/*     KSPSetType(u_ksp, KSPFCG); */
/*     break; */
/*   case CS_PARAM_ITSOL_CG: */
/*     KSPSetType(u_ksp, KSPCG); */
/*     break; */
/*   case CS_PARAM_ITSOL_BICGS:      /\* Improved Bi-CG stab *\/ */
/*     KSPSetType(u_ksp, KSPIBCGS); */
/*     break; */
/*   case CS_PARAM_ITSOL_BICGS2: /\* BiCGstab2 *\/ */
/*     KSPSetType(u_ksp, KSPBCGSL); */
/*     break; */
/*   case CS_PARAM_ITSOL_MUMPS:     /\* Direct solver (factorization) *\/ */
/* #if defined(PETSC_HAVE_MUMPS) */
/*     KSPSetType(u_ksp, KSPPREONLY); */
/*     PCSetType(u_pc, PCLU); */
/*     PCFactorSetMatSolverType(u_pc, MATSOLVERMUMPS); */
/* #else */
/*     bft_error(__FILE__, __LINE__, 0, */
/*               " %s: MUMPS not interfaced with this installation of PETSc.", */
/*               __func__); */
/* #endif */
/*     break; */
/*   case CS_PARAM_ITSOL_GMRES: */
/*     /\* Number of iterations before restarting = 30 (default value)  *\/ */
/*     KSPSetType(u_ksp, KSPGMRES); */
/*     break; */
/*   case CS_PARAM_ITSOL_FGMRES: */
/*     /\* Number of iterations before restarting = 30 (default value)  *\/ */
/*     KSPSetType(u_ksp, KSPFGMRES); */
/*     break; */

/*   default: */
/*     bft_error(__FILE__, __LINE__, 0, "%s: Invalid solver.", __func__); */
/*     break; */

/*   } /\* Switch on solver *\/ */

/*   if (slesp->solver != CS_PARAM_ITSOL_MUMPS) */
/*     PCSetType(u_pc, pc_type); */

/*   /\* Additional settings for the preconditioner *\/ */

/*   switch (slesp->precond) { */

/*   case CS_PARAM_PRECOND_AMG: */
/*     switch (slesp->amg_type) { */

/*     case CS_PARAM_AMG_PETSC_GAMG_V: */
/*       PCGAMGSetType(u_pc, PCGAMGAGG); */
/*       PCGAMGSetNSmooths(u_pc, 1); */
/*       PCMGSetCycleType(u_pc, PC_MG_CYCLE_V); */
/*       _setup_velocity_gamg(); */
/*       break; */

/*     case CS_PARAM_AMG_PETSC_GAMG_W: */
/*       PCGAMGSetType(u_pc, PCGAMGAGG); */
/*       PCGAMGSetNSmooths(u_pc, 1); */
/*       PCMGSetCycleType(u_pc, PC_MG_CYCLE_W); */
/*       _setup_velocity_gamg(); */
/*       break; */

/*     case CS_PARAM_AMG_HYPRE_BOOMER_V: */
/* #if defined(PETSC_HAVE_HYPRE) */
/*       PCHYPRESetType(u_pc, "boomeramg"); */
/*       PetscOptionsSetValue(nullptr, "-pc_hypre_boomeramg_cycle_type","V"); */
/*       _setup_velocity_boomeramg(); */
/* #else */
/*       PCGAMGSetType(u_pc, PCGAMGAGG); */
/*       PCGAMGSetNSmooths(u_pc, 1); */
/*       PCMGSetCycleType(u_pc, PC_MG_CYCLE_V); */
/*       _setup_velocity_gamg(); */
/* #endif */
/*       break; */

/*     case CS_PARAM_AMG_HYPRE_BOOMER_W: */
/* #if defined(PETSC_HAVE_HYPRE) */
/*       PCHYPRESetType(u_pc, "boomeramg"); */
/*       PetscOptionsSetValue(nullptr, "-pc_hypre_boomeramg_cycle_type","W"); */
/*       _setup_velocity_boomeramg(); */
/* #else */
/*       PCGAMGSetType(u_pc, PCGAMGAGG); */
/*       PCGAMGSetNSmooths(u_pc, 1); */
/*       PCMGSetCycleType(u_pc, PC_MG_CYCLE_W); */
/*       _setup_velocity_gamg(); */
/* #endif */
/*       break; */

/*     default: */
/*       bft_error(__FILE__, __LINE__, 0, "%s: Invalid AMG type.", __func__); */
/*       break; */

/*     } /\* AMG type *\/ */
/*     break; */

/*   default: */
/*     break; /\* Nothing else to do *\/ */

/*   } /\* Switch on preconditioner *\/ */

/*   /\* Set tolerance and number of iterations *\/ */

/*   PetscReal _rtol, abstol, dtol; */
/*   PetscInt  _max_it; */
/*   KSPGetTolerances(u_ksp, &_rtol, &abstol, &dtol, &_max_it); */
/*   KSPSetTolerances(u_ksp, */
/*                    rtol,        /\* relative convergence tolerance *\/ */
/*                    abstol,      /\* absolute convergence tolerance *\/ */
/*                    dtol,        /\* divergence tolerance *\/ */
/*                    max_it);     /\* max number of iterations *\/ */

/*   PCSetFromOptions(u_pc); */
/*   PCSetUp(u_pc); */

/*   KSPSetFromOptions(u_ksp); */
/*   KSPSetUp(u_ksp); */
/* } */

#if PETSC_VERSION_GE(3,11,0)
/*----------------------------------------------------------------------------*/
/*!
 * \brief Function pointer: setup hook for setting PETSc solver and
 *        preconditioner.
 *        Case of GKB as a solver.
 *
 * \param[in, out] context     pointer to optional (untyped) value or structure
 * \param[in, out] ksp_struct  pointer to PETSc KSP context
 */
/*----------------------------------------------------------------------------*/

static void
_gkb_hook(void  *context,
          void  *ksp_struct)
{
  /* PETSc structures */

  IS  is1 = nullptr, is2 = nullptr;
  KSP ksp = static_cast<KSP>(ksp_struct);

  cs_fp_exception_disable_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */

  KSPSetType(ksp, KSPPREONLY);

  /* Get the context members */

  cs_saddle_solver_t       *solver = static_cast<cs_saddle_solver_t *>(context);
  const cs_param_saddle_t  *saddlep = solver->param;
  const cs_param_saddle_context_gkb_t *ctxp =
    static_cast<cs_param_saddle_context_gkb_t *>(saddlep->context);

  /* Apply modifications to the KSP structure */

  PC up_pc;

  KSPGetPC(ksp, &up_pc);
  PCSetType(up_pc, PCFIELDSPLIT);
  PCFieldSplitSetType(up_pc, PC_COMPOSITE_GKB);

  PCFieldSplitSetGKBTol(up_pc, saddlep->cvg_param.rtol);
  PCFieldSplitSetGKBMaxit(up_pc, saddlep->cvg_param.n_max_iter);
  PCFieldSplitSetGKBNu(up_pc, ctxp->augmentation_scaling);
  PCFieldSplitSetGKBDelay(up_pc, ctxp->truncation_threshold);

  _build_is_for_fieldsplit(solver, &is1, &is2);

  /* First level block1 | block2 */

  PCFieldSplitSetIS(up_pc, "block1", is1);
  PCFieldSplitSetIS(up_pc, "block2", is2);

  /* Need to call PCSetUp before configuring the second level (Thanks to
     Natacha Bereux) */

  PCSetFromOptions(up_pc);
  PCSetUp(up_pc);
  KSPSetUp(ksp);

  PetscInt  n_split;
  KSP  *up_subksp;
  PCFieldSplitGetSubKSP(up_pc, &n_split, &up_subksp);
  assert(n_split == 2);

  /* Set the (1,1)-block
   * One needs a "non const" version of the structure for the setup */

  cs_equation_param_t  *eqp =
    cs_equation_param_by_name(saddlep->block11_sles_param->name);
  assert(eqp != nullptr);

  cs_param_sles_t  *block11_slesp = eqp->sles_param;

  cs_param_sles_setup_petsc_ksp("block1", block11_slesp, up_subksp[0]);

  /* Apply modifications to the KSP structure */

  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);

  PetscFree(up_subksp);
  ISDestroy(&is2);
  ISDestroy(&is1);

  cs_fp_exception_restore_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */
}
#endif  /* PETSC VERSION >= 3.11 */
#endif  /* HAVE_PETSC */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function pointer: setup hook for setting PETSc
 *        Case of Notay's algebraic transformation.
 *        This relies on the following article:
 *        "Algebraic multigrid for Stokes equations", Y. Notay (2017)
 *        SIAM J. Sci. Comput., Vol. 39 (5), pp 88-111
 *
 * \param[in, out] context    pointer to optional (untyped) value or structure
 * \param[in, out] ksp_struct pointer to PETSc KSP context
 */
/*----------------------------------------------------------------------------*/

static void
_notay_hook(void  *context,
            void  *ksp_struct)
{
  cs_saddle_solver_t *solver = static_cast<cs_saddle_solver_t *>(context);

  const cs_param_saddle_t  *saddlep = solver->param;
  const cs_param_saddle_context_notay_t *ctxp =
    static_cast<const cs_param_saddle_context_notay_t *>(saddlep->context);
  cs_param_sles_t *b11_slesp =
    const_cast<cs_param_sles_t *>(saddlep->block11_sles_param);

#if defined(HAVE_PETSC)
  KSP ksp = static_cast<KSP>(ksp_struct);

  cs_fp_exception_disable_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */

  if (cs_glob_n_ranks > 1) {
    cs_base_warn(__FILE__, __LINE__);
    cs_log_printf(CS_LOG_DEFAULT,
                  " %s (Eq. %s) Warning: Algo. not tested in parallel.\n",
                  __func__, b11_slesp->name);
  }

  /* Build IndexSet structures to extract block matrices */

  IS  is1 = nullptr, is2 = nullptr;
  _build_is_for_fieldsplit(solver, &is1, &is2);

  Mat Amat, Amat_nest;
  KSPGetOperators(ksp, &Amat, nullptr);

  /* Retrieve blocks */

  Mat A11, A12, A21;
  MatCreateSubMatrix(Amat, is1, is1, MAT_INITIAL_MATRIX, &A11);
  MatCreateSubMatrix(Amat, is1, is2, MAT_INITIAL_MATRIX, &A12);
  MatCreateSubMatrix(Amat, is2, is1, MAT_INITIAL_MATRIX, &A21);

  PetscInt n1, n2;
  MatGetSize(A12, &n1, &n2);

  /* Define diag = inv(diag(A11)) */

  Vec D11;
  VecCreate(PETSC_COMM_WORLD, &D11);
  VecSetSizes(D11, PETSC_DECIDE, n1);
  VecSetType(D11, VECMPI);
  MatGetDiagonal(A11, D11);
  VecReciprocal(D11);

  const PetscReal  alpha = ctxp->scaling_coef;
  if (fabs(alpha - 1.0) > 0)
    VecScale(D11, alpha);

  /* Computing new blocks for the transformed system */

  PetscScalar one = 1.0;

  /* Compute the block (1, 2)
   * First step: T11 <- Id - A11*inv(diag(A11)) */

  Mat T11;
  MatConvert(A11, MATSAME, MAT_INITIAL_MATRIX, &T11);

  MatDiagonalScale(T11, nullptr, D11); /* left scaling = nullptr;
                                       right scaling = D11 */
  MatScale(T11, -one);

  Vec ones;
  VecCreate(PETSC_COMM_WORLD, &ones);
  VecSetSizes(ones, PETSC_DECIDE, n1);
  VecSetType(ones, VECMPI);
  VecSet(ones, one);
  MatDiagonalSet(T11, ones, ADD_VALUES);

  /* Second step: T12 = T11*A12 = (Id - A11*inv(D_A11))*A12 */

  Mat T12;
  MatMatMult(T11, A12, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &T12);

  /* Compute the A22 block
   * A22 = A21*inv(diag(A11))*A12 */

  Mat T21, A22;
  MatConvert(A21, MATSAME, MAT_INITIAL_MATRIX, &T21);
  MatDiagonalScale(T21, nullptr, D11); /* T21 <- A21*inv(diag(A11)) */
  MatMatMult(T21, A12, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &A22);

  /* Partial free */

  VecDestroy(&D11);
  VecDestroy(&ones);
  MatDestroy(&A12);
  MatDestroy(&T11);
  MatDestroy(&T21);

  /* Compute the (2, 1)-block
   *   A21 <- -1.0*A21 */

  MatScale(A21, -1.0);

  /* Update blocks and assemble Amat */

  Mat subA[4] = {A11, T12, A21, A22};

  MatCreateNest(PETSC_COMM_WORLD, 2, nullptr, 2, nullptr, subA, &Amat_nest);
  MatConvert(Amat_nest, MATMPIAIJ, MAT_INITIAL_MATRIX, &Amat);

  KSPSetOperators(ksp, Amat, Amat);

  /* Partial free */

  MatDestroy(&Amat_nest);
  MatDestroy(&A11);
  MatDestroy(&A21);
  MatDestroy(&A22);
  MatDestroy(&T12);
  ISDestroy(&is1);
  ISDestroy(&is2);
  MatDestroy(&Amat);

  PC  saddle_pc;
  KSPGetPC(ksp, &saddle_pc);

  /* Set the main solver and main options for the preconditioner */

  cs_param_sles_setup_petsc_ksp(saddlep->name, b11_slesp, ksp);

  cs_fp_exception_restore_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */
#else
  bft_error(__FILE__, __LINE__, 0,
            "%s: PETSc is required for solving \"%s\"\n"
            " Please modify your settings/build code_saturne with PETSc.",
            __func__, saddlep->name);
#endif  /* HAVE_PETSC */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Do the setup stage in case of GKB algorithm
 *
 * \param[in, out] solver         pointer to a saddle-point solver structure
 * \param[in, out] saddlep        set of parameters for solving a saddle-point
 * \param[in, out] block11_slesp  set of parameters for the (1,1)-block system
 */
/*----------------------------------------------------------------------------*/

static void
_gkb_setup(cs_saddle_solver_t  *solver,
           cs_param_saddle_t   *saddlep,
           cs_param_sles_t     *block11_slesp)
{
  assert(saddlep->solver == CS_PARAM_SADDLE_SOLVER_GKB);

  switch (saddlep->solver_class) {

  case CS_PARAM_SOLVER_CLASS_PETSC:
    {
      /* Golub-Kahan Bi-diagonalization is available starting from the 3.11
         version of PETSc */
#if defined(HAVE_PETSC)
#if PETSC_VERSION_GE(3,11,0)
      cs_sles_petsc_init();

      cs_sles_petsc_define(block11_slesp->field_id,
                           nullptr,
                           MATMPIAIJ,
                           _gkb_hook,
                           (void *)solver); /* context structure for PETSc */
#else
      bft_error(__FILE__, __LINE__, 0,
                "%s: Invalid GKB solver in PETSc\n"
                " PETSc 3.11.x or greater is required.\n", __func__);
#endif
#endif  /* HAVE PETSC */
    }
    break;

  default: /* in-house version of the GKB */
    {
      int  ierr = cs_param_sles_setup(true, block11_slesp);

      if (ierr < 0)
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Error %d detected during the setup of the (1,1)-block"
                  " of the \"%s\" saddle-point problem with GKB.\n"
                  " Please modify your settings.",
                  __func__, ierr, saddlep->name);

      /* SLES for the initial transformation of the right hand-side */

      cs_param_saddle_context_gkb_t *ctxp =
        static_cast<cs_param_saddle_context_gkb_t *>(saddlep->context);

      if (ctxp->dedicated_init_sles) {

        ierr = cs_param_sles_setup(false, ctxp->init_sles_param);

        if (ierr < 0)
          bft_error(__FILE__, __LINE__, 0,
                    "%s: Error %d detected during the setup of the initial SLES"
                    " of a saddle-point problem with GKB.\n"
                    " Please modify your settings.",
                    __func__, ierr);

      }
    }
    break;

  } /* Class of solver to consider */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Do the setup stage when the Notay's algebraic transformation is
 *        chosen
 *
 * \param[in, out] solver         pointer to a saddle-point solver structure
 * \param[in, out] saddlep        set of parameters for solving a saddle-point
 * \param[in, out] block11_slesp  set of parameters for the (1,1)-block system
 */
/*----------------------------------------------------------------------------*/

static void
_notay_transformation_setup(cs_saddle_solver_t  *solver,
                            cs_param_saddle_t   *saddlep,
                            cs_param_sles_t     *block11_slesp)
{
  assert(saddlep->solver == CS_PARAM_SADDLE_SOLVER_NOTAY_TRANSFORM);

  switch (saddlep->solver_class) {

  case CS_PARAM_SOLVER_CLASS_PETSC:
#if defined(HAVE_PETSC)
    cs_sles_petsc_init();

    cs_sles_petsc_define(block11_slesp->field_id,
                         nullptr,
                         MATMPIAIJ,
                         _notay_hook,
                         (void *)solver); /* context structure for PETSc */
#else
    bft_error(__FILE__, __LINE__, 0,
              "%s: PETSc is required for solving \"%s\"\n",
              " Please modify your settings/build code_saturne with PETSc.",
              __func__, saddlep->name);
#endif  /* HAVE PETSC */
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid solver class.\n"
              " PETSc is requested for solving \"%s\"\n"
              " Please modify your settings.",
              __func__, saddlep->name);
    break;

  } /* Class of solver to consider */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the SLES structure(s) related to the resolution or the
 *        construction of the Schur complement approximation
 *
 * \param[in, out] saddlep  set of parameters for solving a saddle-point sys.
 *
 * \return an error code (< 0 if a problem is encountered, 0 otherwise)
 */
/*----------------------------------------------------------------------------*/

static int
_schur_complement_setup(cs_param_saddle_t  *saddlep)
{
  int ierr = 0; /* OK */

  if (saddlep->schur_approx == CS_PARAM_SADDLE_SCHUR_NONE     ||
      saddlep->schur_approx == CS_PARAM_SADDLE_SCHUR_IDENTITY ||
      saddlep->schur_approx == CS_PARAM_SADDLE_SCHUR_MASS_SCALED)
    return ierr; /* Nothing to do */

  cs_param_sles_t  *schur_slesp = saddlep->schur_sles_param;
  assert(schur_slesp != nullptr);

  ierr = cs_param_sles_setup(false, schur_slesp);

  if (ierr < 0) {
    cs_base_warn(__FILE__, __LINE__);
    cs_log_printf(CS_LOG_WARNINGS,
                  "%s: Problem detected when settings the system \"%s\""
                  " related to a Schur complement\n",
                  __func__, schur_slesp->name);
    return ierr;
  }

  switch (saddlep->schur_approx) {

  case CS_PARAM_SADDLE_SCHUR_LUMPED_INVERSE:
  case CS_PARAM_SADDLE_SCHUR_MASS_SCALED_LUMPED_INVERSE:
    {
      cs_param_sles_t  *xtra_slesp =
        cs_param_saddle_get_xtra_sles_param(saddlep);
      assert(xtra_slesp != nullptr);

      ierr = cs_param_sles_setup(false, xtra_slesp);

      if (ierr < 0) {

        cs_base_warn(__FILE__, __LINE__);
        cs_log_printf(CS_LOG_WARNINGS,
                      "%s: Problem detected when settings the system \"%s\""
                      " related to the Schur complement (extra system)\n",
                      __func__, xtra_slesp->name);

        return ierr;
      }
    }
    break;

  default:
    /* Do nothing */
    break;

  } /* Switch on the type of approximation for the Schur complement */

  return ierr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define at least a cs_sles_t structure in accordance with the settings
 *        of a cs_param_saddle_t structure (SLES = Sparse Linear Equation
 *        Solver). Some \ref cs_param_sles_t structure are embedded in the
 *        strategy to solve a saddle-point problem.
 *
 * \param[in, out] solver         pointer to a saddle-point solver structure
 * \param[in, out] saddlep        set of parameters for solving a saddle-point
 * \param[in, out] block11_slesp  set of parameters for the (1,1)-block system
 */
/*----------------------------------------------------------------------------*/

static void
_setup(cs_saddle_solver_t  *solver,
       cs_param_saddle_t   *saddlep,
       cs_param_sles_t     *block11_slesp)
{
  int ierr;

  /* Setup according to the type of solver */

  switch (saddlep->solver) {

  case CS_PARAM_SADDLE_SOLVER_ALU:
    /* -------------------------- */
    {
      cs_param_saddle_context_alu_t *ctxp =
        static_cast<cs_param_saddle_context_alu_t *>(saddlep->context);

      ierr = cs_param_sles_setup(true, block11_slesp);

      if (ierr < 0)
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Error %d detected during the setup of the (1,1)-block"
                  " of a saddle-point problem with ALU.\n"
                  " Please modify your settings.",
                  __func__, ierr);

      /* SLES for the initial transformation of the right hand-side */

      if (ctxp->dedicated_init_sles) {

        ierr = cs_param_sles_setup(false, ctxp->init_sles_param);

        if (ierr < 0)
          bft_error(__FILE__, __LINE__, 0,
                    "%s: Error %d detected during the setup of the initial SLES"
                    " of a saddle-point problem with ALU.\n"
                    " Please modify your settings.",
                    __func__, ierr);

      }
    }
    break;

  case CS_PARAM_SADDLE_SOLVER_FGMRES:
    /* ----------------------------- */
    //_petsc_fmgres_setup(solver, saddlep, block11_slesp);
    break;

  case CS_PARAM_SADDLE_SOLVER_GCR:
  case CS_PARAM_SADDLE_SOLVER_MINRES:
    /* ----------------------------- */
    ierr = cs_param_sles_setup(true, block11_slesp);

    if (ierr < 0)
      bft_error(__FILE__, __LINE__, 0,
                "%s: Error %d detected during the setup of the (1,1)-block"
                " of the \"%s\" saddle-point problem with GCR/MINRES.\n"
                " Please check your settings.",
                __func__, ierr, saddlep->name);

    /* Handle additional SLES systems */

    ierr = _schur_complement_setup(saddlep);

    if (ierr < 0)
      bft_error(__FILE__, __LINE__, 0,
                "%s: Error %d detected during the setup related to the Schur"
                " complement associated to the \"%s\" saddle-point problem.\n"
                " Please check your settings.",
                __func__, ierr, saddlep->name);
    break;

  case CS_PARAM_SADDLE_SOLVER_MUMPS:
    /* ---------------------------- */
    if (block11_slesp->solver != CS_PARAM_SOLVER_MUMPS) {

      block11_slesp->solver = CS_PARAM_SOLVER_MUMPS;
      block11_slesp->precond = CS_PARAM_PRECOND_NONE;
      block11_slesp->solver_class = CS_PARAM_SOLVER_CLASS_MUMPS;

    }

    ierr = cs_param_sles_setup(true, block11_slesp);
    if (ierr < 0)
      bft_error(__FILE__, __LINE__, 0,
                "%s: Error %d detected during the setup"
                " of the \"%s\" saddle-point problem with MUMPS.\n"
                " Please check your settings.",
                __func__, ierr, saddlep->name);
    break; /* MUMPS */

  case CS_PARAM_SADDLE_SOLVER_GKB:
    /* -------------------------- */
    _gkb_setup(solver, saddlep, block11_slesp);

    break; /* GKB */

  case CS_PARAM_SADDLE_SOLVER_NOTAY_TRANSFORM:
    /* --------------------------------------- */
    _notay_transformation_setup(solver, saddlep, block11_slesp);
    break;

  case CS_PARAM_SADDLE_SOLVER_UZAWA_CG:
    /* ------------------------------- */
    ierr = cs_param_sles_setup(true, block11_slesp);

    if (ierr < 0)
      bft_error(__FILE__, __LINE__, 0,
                "%s: Error %d detected during the setup of the (1,1)-block"
                " of the \"%s\" saddle-point problem with Uzawa-CG.\n"
                " Please check your settings.",
                __func__, ierr, saddlep->name);

    /* Handle additional SLES systems */

    ierr = _schur_complement_setup(saddlep);

    if (ierr < 0)
      bft_error(__FILE__, __LINE__, 0,
                "%s: Error %d detected during the setup related to the Schur"
                " complement associated to the \"%s\" saddle-point problem.\n"
                " Please check your settings.",
                __func__, ierr, saddlep->name);
    break;

  case CS_PARAM_SADDLE_SOLVER_SIMPLE:
    /* ------------------------------- */
    ierr = cs_param_sles_setup(true, block11_slesp);

    if (ierr < 0)
      bft_error(__FILE__, __LINE__, 0,
                "%s: Error %d detected during the setup of the (1,1)-block"
                " of the \"%s\" saddle-point problem with SIMPLE.\n"
                " Please check your settings.",
                __func__, ierr, saddlep->name);

    /* Handle additional SLES systems */

    ierr = _schur_complement_setup(saddlep);

    if (ierr < 0)
      bft_error(__FILE__, __LINE__, 0,
                "%s: Error %d detected during the setup related to the Schur"
                " complement associated to the \"%s\" saddle-point problem.\n"
                " Please check your settings.",
                __func__, ierr, saddlep->name);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Unknow type of saddle-point solver.\n"
              " Please check your settings for the \"%s\" saddle-point system.",
              __func__, saddlep->name);
    break;

  } /* End of switch on class of solvers */
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the SLES structures (potentially several) related to a
 *        saddle-point system
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_setup_sles(void)
{
  int n = cs_saddle_solver_get_n_systems();

  for (int id = 0; id < n; id++) {

    cs_saddle_solver_t  *solver = cs_saddle_solver_by_id(id);
    assert(solver != nullptr);
    const cs_param_saddle_t  *const_saddlep = solver->param;

    if (const_saddlep->solver == CS_PARAM_SADDLE_SOLVER_NONE)
      continue;

    /* One needs a "non const" version of the structures for the setup */

    cs_equation_param_t  *eqp =
      cs_equation_param_by_name(const_saddlep->block11_sles_param->name);
    assert(eqp != nullptr);

    cs_param_saddle_t  *saddlep = eqp->saddle_param;
    cs_param_sles_t  *block11_slesp = eqp->sles_param;

    _setup(solver, saddlep, block11_slesp);

  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
