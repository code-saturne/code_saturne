/*============================================================================
 * Routines to handle the SLES settings
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

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#if defined(HAVE_PETSC)
#include <petsc.h>
#include <petscconf.h> /* Useful to know if HYPRE is accessible through PETSc */
#include <petscversion.h>
#endif

#if defined(HAVE_HYPRE)
#include <HYPRE_krylov.h>
#include <HYPRE_parcsr_ls.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_error.h>
#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_fp_exception.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_multigrid.h"
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

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return true if the prescribed solver implies a symmetric linear
 *        system
 */
/*----------------------------------------------------------------------------*/

static inline bool
_system_should_be_sym(cs_param_itsol_type_t   solver)
{
  switch (solver) {

  case CS_PARAM_ITSOL_CG:
  case CS_PARAM_ITSOL_FCG:
  case CS_PARAM_ITSOL_GKB_CG:
  case CS_PARAM_ITSOL_GKB_GMRES:
  case CS_PARAM_ITSOL_MINRES:
    return true;

  default:
    return false; /* Assume that the system is not symmetric */
  }
}

#if defined(HAVE_PETSC)
/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the command line option for PETSc
 *
 * \param[in]      use_prefix    need a prefix
 * \param[in]      prefix        optional prefix
 * \param[in]      keyword       command keyword
 * \param[in]      keyval        command value
 */
/*----------------------------------------------------------------------------*/

static inline void
_petsc_cmd(bool          use_prefix,
           const char   *prefix,
           const char   *keyword,
           const char   *keyval)
{
  char  cmd_line[128];

  if (use_prefix)
    sprintf(cmd_line, "-%s_%s", prefix, keyword);
  else
    sprintf(cmd_line, "-%s", keyword);

#if PETSC_VERSION_GE(3,7,0)
  PetscOptionsSetValue(NULL, cmd_line, keyval);
#else
  PetscOptionsSetValue(cmd_line, keyval);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Predefined settings for a block ILU(0) with PETSc
 *
 * \param[in]      prefix        prefix name associated to the current SLES
 */
/*----------------------------------------------------------------------------*/

static inline void
_petsc_bilu0_hook(const char              *prefix)
{
  assert(prefix != NULL);

  _petsc_cmd(true, prefix, "pc_type", "bjacobi");
  _petsc_cmd(true, prefix, "pc_jacobi_blocks", "1");
  _petsc_cmd(true, prefix, "sub_ksp_type", "preonly");
  _petsc_cmd(true, prefix, "sub_pc_type", "ilu");
  _petsc_cmd(true, prefix, "sub_pc_factor_level", "0");
  _petsc_cmd(true, prefix, "sub_pc_factor_reuse_ordering", "");

  /* If one wants to optimize the memory consumption */
  /* _petsc_cmd(true, prefix, "sub_pc_factor_in_place", ""); */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Predefined settings for a block ICC(0) with PETSc
 *
 * \param[in]      prefix        prefix name associated to the current SLES
 */
/*----------------------------------------------------------------------------*/

static inline void
_petsc_bicc0_hook(const char              *prefix)
{
  assert(prefix != NULL);

  _petsc_cmd(true, prefix, "pc_type", "bjacobi");
  _petsc_cmd(true, prefix, "pc_jacobi_blocks", "1");
  _petsc_cmd(true, prefix, "sub_ksp_type", "preonly");
  _petsc_cmd(true, prefix, "sub_pc_type", "icc");
  _petsc_cmd(true, prefix, "sub_pc_factor_level", "0");
  _petsc_cmd(true, prefix, "sub_pc_factor_reuse_ordering", "");

  /* If one wants to optimize the memory consumption */
  /* _petsc_cmd(true, prefix, "sub_pc_factor_in_place", ""); */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Predefined settings for a block SSOR with PETSc
 *
 * \param[in]      prefix        prefix name associated to the current SLES
 */
/*----------------------------------------------------------------------------*/

static inline void
_petsc_bssor_hook(const char              *prefix)
{
  assert(prefix != NULL);

  _petsc_cmd(true, prefix, "pc_type", "bjacobi");
  _petsc_cmd(true, prefix, "pc_jacobi_blocks", "1");
  _petsc_cmd(true, prefix, "sub_ksp_type", "preonly");
  _petsc_cmd(true, prefix, "sub_pc_type", "sor");
  _petsc_cmd(true, prefix, "sub_pc_sor_symmetric", "");
  _petsc_cmd(true, prefix, "sub_pc_sor_local_symmetric", "");
  _petsc_cmd(true, prefix, "sub_pc_sor_omega", "1.5");
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Predefined settings for GAMG as a preconditioner even if another
 *        settings have been defined. One assumes that one really wants to use
 *        GAMG (may be HYPRE is not available)
 *
 * \param[in]      prefix        prefix name associated to the current SLES
 * \param[in]      slesp         pointer to a set of SLES parameters
 * \param[in]      is_symm       the linear system to solve is symmetric
 * \param[in, out] pc            pointer to a PETSc preconditioner
 */
/*----------------------------------------------------------------------------*/

static inline void
_petsc_pcgamg_hook(const char              *prefix,
                   const cs_param_sles_t   *slesp,
                   bool                     is_symm,
                   PC                       pc)
{
  assert(prefix != NULL);
  assert(slesp != NULL);
  assert(slesp->precond == CS_PARAM_PRECOND_AMG);

  /* Remark: -pc_gamg_reuse_interpolation
   *
   * Reuse prolongation when rebuilding algebraic multigrid
   * preconditioner. This may negatively affect the convergence rate of the
   * method on new matrices if the matrix entries change a great deal, but
   * allows rebuilding the preconditioner quicker. (default=false)
   */

  _petsc_cmd(true, prefix, "pc_gamg_reuse_interpolation", "true");

  /* Remark: -pc_gamg_sym_graph
   * Symmetrize the graph before computing the aggregation. Some algorithms
   * require the graph be symmetric (default=false)
   */

  _petsc_cmd(true, prefix, "pc_gamg_sym_graph", "true");

  /* Set smoothers (general settings, i.e. not depending on the symmetry or not
     of the linear system to solve) */

  _petsc_cmd(true, prefix, "mg_levels_ksp_type", "richardson");
  _petsc_cmd(true, prefix, "mg_levels_ksp_max_it", "1");
  _petsc_cmd(true, prefix, "mg_levels_ksp_norm_type", "none");
  _petsc_cmd(true, prefix, "mg_levels_ksp_richardson_scale", "1.0");

  /* Do not build a coarser level if one reaches the following limit */

  _petsc_cmd(true, prefix, "pc_gamg_coarse_eq_limit", "100");

  /* In parallel computing, migrate data to another rank if the grid has less
     than 200 rows */

  if (cs_glob_n_ranks > 1) {

    _petsc_cmd(true, prefix, "pc_gamg_repartition", "true");
    _petsc_cmd(true, prefix, "pc_gamg_process_eq_limit", "200");

  }
  else {

    _petsc_cmd(true, prefix, "mg_coarse_ksp_type", "preonly");
    _petsc_cmd(true, prefix, "mg_coarse_pc_type", "tfs");

  }

  /* Settings depending on the symmetry or not of the linear system to solve */

  if (is_symm) {

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

    _petsc_cmd(true, prefix, "pc_gamg_agg_nsmooths", "2");
    _petsc_cmd(true, prefix, "pc_gamg_square_graph", "2");
    _petsc_cmd(true, prefix, "pc_gamg_threshold", "0.08");

    if (cs_glob_n_ranks > 1) {

      _petsc_cmd(true, prefix, "mg_levels_pc_type", "bjacobi");
      _petsc_cmd(true, prefix, "mg_levels_pc_jacobi_blocks", "1");
      _petsc_cmd(true, prefix, "mg_levels_sub_ksp_type", "preonly");
      _petsc_cmd(true, prefix, "mg_levels_sub_pc_type", "sor");
      _petsc_cmd(true, prefix, "mg_levels_sub_pc_sor_local_symmetric", "");
      _petsc_cmd(true, prefix, "mg_levels_sub_pc_sor_omega", "1.5");

    }
    else { /* serial run */

      _petsc_cmd(true, prefix, "mg_levels_pc_type", "sor");
      _petsc_cmd(true, prefix, "mg_levels_pc_sor_local_symmetric", "");
      _petsc_cmd(true, prefix, "mg_levels_pc_sor_omega", "1.5");

    }

  }
  else { /* Not a symmetric linear system */

    /* Number of smoothing steps to use with smooth aggregation (default=1) */

    _petsc_cmd(true, prefix, "pc_gamg_agg_nsmooths", "0");
    _petsc_cmd(true, prefix, "pc_gamg_square_graph", "0");
    _petsc_cmd(true, prefix, "pc_gamg_threshold", "0.06");

    _petsc_cmd(true, prefix, "mg_levels_pc_type", "bjacobi");
    _petsc_cmd(true, prefix, "mg_levels_pc_bjacobi_blocks", "1");
    _petsc_cmd(true, prefix, "mg_levels_sub_ksp_type", "preonly");
    _petsc_cmd(true, prefix, "mg_levels_sub_pc_type", "ilu");
    _petsc_cmd(true, prefix, "mg_levels_sub_pc_factor_levels", "0");

  }

  /* After command line options, switch to PETSc setup functions */

  PCSetType(pc, PCGAMG);
  PCGAMGSetType(pc, PCGAMGAGG);
  PCGAMGSetNSmooths(pc, 1);
  PCSetUp(pc);

  switch (slesp->amg_type) {

  case CS_PARAM_AMG_PETSC_GAMG_V:
  case CS_PARAM_AMG_PETSC_PCMG:
  case CS_PARAM_AMG_HYPRE_BOOMER_V:
    PCMGSetCycleType(pc, PC_MG_CYCLE_V);
    break;

  case CS_PARAM_AMG_PETSC_GAMG_W:
  case CS_PARAM_AMG_HYPRE_BOOMER_W:
    PCMGSetCycleType(pc, PC_MG_CYCLE_W);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid type of AMG for SLES %s\n",
              __func__, slesp->name);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Predefined settings for BoomerAMG in HYPRE as a preconditioner
 *
 * \param[in]      prefix        prefix name associated to the current SLES
 * \param[in]      slesp         pointer to a set of SLES parameters
 * \param[in]      is_symm       the linear system to solve is symmetric
 * \param[in, out] pc            pointer to a PETSc preconditioner
 */
/*----------------------------------------------------------------------------*/

static inline void
_petsc_pchypre_hook(const char              *prefix,
                    const cs_param_sles_t   *slesp,
                    bool                     is_symm,
                    PC                       pc)
{
  CS_UNUSED(is_symm);

  assert(prefix != NULL);
  assert(slesp != NULL);
  assert(slesp->precond == CS_PARAM_PRECOND_AMG);

  PCSetType(pc, PCHYPRE);
  PCHYPRESetType(pc, "boomeramg");

  switch (slesp->amg_type) {

  case CS_PARAM_AMG_HYPRE_BOOMER_V:
    _petsc_cmd(true, prefix, "pc_hypre_boomeramg_cycle_type", "V");
    break;

  case CS_PARAM_AMG_HYPRE_BOOMER_W:
    _petsc_cmd(true, prefix, "pc_hypre_boomeramg_cycle_type", "W");
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid type of AMG for SLES %s\n",
              __func__, slesp->name);
  }

  /* From HYPRE documentation: https://hypre.readthedocs.io/en/lastest
   *
   * for three-dimensional diffusion problems, it is recommended to choose a
   * lower complexity coarsening like HMIS or PMIS (coarsening 10 or 8) and
   * combine it with a distance-two interpolation (interpolation 6 or 7), that
   * is also truncated to 4 or 5 elements per row. Additional reduction in
   * complexity and increased scalability can often be achieved using one or
   * two levels of aggressive coarsening.
   */

  /* _petsc_cmd(true, prefix, "pc_hypre_boomeramg_grid_sweeps_down","2"); */
  /* _petsc_cmd(true, prefix, "pc_hypre_boomeramg_grid_sweeps_up","2"); */
  /* _petsc_cmd(true, prefix, "pc_hypre_boomeramg_smooth_type","Euclid"); */

  /* Remark: fcf-jacobi or l1scaled-jacobi (or chebyshev) as up/down smoothers
     can be a good choice
  */

  _petsc_cmd(true, prefix,
             // "pc_hypre_boomeramg_relax_type_down", "l1scaled-SOR/Jacobi");
             "pc_hypre_boomeramg_relax_type_down", "symmetric-SOR/Jacobi");

  _petsc_cmd(true, prefix,
             // "pc_hypre_boomeramg_relax_type_up", "l1scaled-SOR/Jacobi");
             "pc_hypre_boomeramg_relax_type_up", "symmetric-SOR/Jacobi");

  if (slesp->solver == CS_PARAM_ITSOL_CG ||
      slesp->solver == CS_PARAM_ITSOL_FCG ||
      slesp->solver == CS_PARAM_ITSOL_MINRES)
    _petsc_cmd(true, prefix,
               "pc_hypre_boomeramg_relax_type_coarse", "CG");
  else
    _petsc_cmd(true, prefix,
               "pc_hypre_boomeramg_relax_type_coarse", "Gaussian-elimination");

  /* Note that the default coarsening is HMIS in HYPRE */

  _petsc_cmd(true, prefix, "pc_hypre_boomeramg_coarsen_type", "HMIS");

  /* Note that the default interpolation is extended+i interpolation truncated
   * to 4 elements per row. Using 0 means there is no limitation.
   * good choices are: ext+i-cc, ext+i, FF1
   */

  _petsc_cmd(true, prefix, "pc_hypre_boomeramg_interp_type", "ext+i-cc");
  _petsc_cmd(true, prefix, "pc_hypre_boomeramg_P_max","8");

  /* Number of levels (starting from the finest one) on which one applies an
     aggressive coarsening */

  _petsc_cmd(true, prefix, "pc_hypre_boomeramg_agg_nl","2");

  /* Number of paths for aggressive coarsening (default = 1) */

  _petsc_cmd(true, prefix, "pc_hypre_boomeramg_agg_num_paths","2");

  /* For best performance, it might be necessary to set certain parameters,
   * which will affect both coarsening and interpolation. One important
   * parameter is the strong threshold.  The default value is 0.25, which
   * appears to be a good choice for 2-dimensional problems and the low
   * complexity coarsening algorithms. For 3-dimensional problems a better
   * choice appears to be 0.5, when using the default coarsening
   * algorithm. However, the choice of the strength threshold is problem
   * dependent.
   */

  _petsc_cmd(true, prefix, "pc_hypre_boomeramg_strong_threshold","0.5");
  _petsc_cmd(true, prefix, "pc_hypre_boomeramg_no_CF","");
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set command line options for PC according to the kind of
 *        preconditionner
 *
 * \param[in, out] slesp    set of parameters for the linear algebra
 * \param[in, out] ksp      PETSc solver structure
 */
/*----------------------------------------------------------------------------*/

static void
_petsc_set_pc_type(cs_param_sles_t   *slesp,
                   KSP                ksp)
{
  if (cs_param_sles_is_mumps_set(slesp->solver))
    return; /* Direct solver: Nothing to do at this stage */

  PC  pc;
  KSPGetPC(ksp, &pc);

  switch (slesp->precond) {

  case CS_PARAM_PRECOND_NONE:
    PCSetType(pc, PCNONE);
    break;

  case CS_PARAM_PRECOND_DIAG:
    PCSetType(pc, PCJACOBI);
    break;

  case CS_PARAM_PRECOND_BJACOB_ILU0:
    if (slesp->solver_class == CS_PARAM_SLES_CLASS_HYPRE) {
#if defined(PETSC_HAVE_HYPRE)
      PCSetType(pc, PCHYPRE);
      PCHYPRESetType(pc, "euclid");
      _petsc_cmd(true, slesp->name, "pc_euclid_level", "0");
#else
      _petsc_bilu0_hook(slesp->name);
#endif
    }
    else
      _petsc_bilu0_hook(slesp->name);
    break;

  case CS_PARAM_PRECOND_BJACOB_SGS:
    _petsc_bssor_hook(slesp->name);
    break;

  case CS_PARAM_PRECOND_SSOR:
    if (cs_glob_n_ranks > 1) {  /* Switch to a block version */

      slesp->precond = CS_PARAM_PRECOND_BJACOB_SGS;
      cs_base_warn(__FILE__, __LINE__);
      cs_log_printf(CS_LOG_DEFAULT,
                    " %s: System %s: Modify the requested preconditioner to"
                    " enable a parallel computation with PETSC.\n"
                    " Switch to a block jacobi preconditioner.\n",
                    __func__, slesp->name);

      _petsc_bssor_hook(slesp->name);

    }
    else { /* Serial computation */
      PCSetType(pc, PCSOR);
      PCSORSetSymmetric(pc, SOR_SYMMETRIC_SWEEP);
    }
    break;

  case CS_PARAM_PRECOND_ICC0:
    if (cs_glob_n_ranks > 1) { /* Switch to a block version */

      cs_base_warn(__FILE__, __LINE__);
      cs_log_printf(CS_LOG_DEFAULT,
                    " %s: System %s: Modify the requested preconditioner to"
                    " enable a parallel computation with PETSC.\n"
                    " Switch to a block jacobi preconditioner.\n",
                    __func__, slesp->name);

      _petsc_bicc0_hook(slesp->name);

    }
    else {
      PCSetType(pc, PCICC);
      PCFactorSetLevels(pc, 0);
    }
    break;

  case CS_PARAM_PRECOND_ILU0:
    if (slesp->solver_class == CS_PARAM_SLES_CLASS_HYPRE) {
#if defined(PETSC_HAVE_HYPRE)
      /* Euclid is a parallel version of the ILU(0) factorisation */
      PCSetType(pc, PCHYPRE);
      PCHYPRESetType(pc, "euclid");
      _petsc_cmd(true, slesp->name, "pc_euclid_level", "0");
#else
      _petsc_bilu0_hook(slesp->name);
      if (cs_glob_n_ranks > 1)  /* Switch to a block version */
        slesp->precond = CS_PARAM_PRECOND_BJACOB_ILU0;
#endif
    }
    else {
      _petsc_bilu0_hook(slesp->name);
      if (cs_glob_n_ranks > 1) { /* Switch to a block version */

        slesp->precond = CS_PARAM_PRECOND_BJACOB_ILU0;
        cs_base_warn(__FILE__, __LINE__);
        cs_log_printf(CS_LOG_DEFAULT,
                      " %s: System %s: Modify the requested preconditioner to"
                      " enable a parallel computation with PETSC.\n"
                      " Switch to a block jacobi preconditioner.\n",
                      __func__, slesp->name);
      }
    }
    break;

  case CS_PARAM_PRECOND_LU:
#if defined(PETSC_HAVE_MUMPS)
  case CS_PARAM_PRECOND_MUMPS:
    _petsc_cmd(true, slesp->name, "pc_type", "lu");
    _petsc_cmd(true, slesp->name, "pc_factor_mat_solver_type", "mumps");
#else
    if (cs_glob_n_ranks == 1)
      _petsc_cmd(true, slesp->name, "pc_type", "lu");
    else { /* Switch to a block version (sequential in each block) */
      _petsc_cmd(true, slesp->name, "pc_type", "bjacobi");
      _petsc_cmd(true, slesp->name, "pc_jacobi_blocks", "1");
      _petsc_cmd(true, slesp->name, "sub_ksp_type", "preonly");
      _petsc_cmd(true, slesp->name, "sub_pc_type", "lu");
    }
#endif
    break;

  case CS_PARAM_PRECOND_AMG:
    {
      bool  is_symm = _system_should_be_sym(slesp->solver);

      switch (slesp->amg_type) {

      case CS_PARAM_AMG_PETSC_GAMG_V:
      case CS_PARAM_AMG_PETSC_GAMG_W:
      case CS_PARAM_AMG_PETSC_PCMG:
        _petsc_pcgamg_hook(slesp->name, slesp, is_symm, pc);
        break;

      case CS_PARAM_AMG_HYPRE_BOOMER_V:
      case CS_PARAM_AMG_HYPRE_BOOMER_W:
#if defined(PETSC_HAVE_HYPRE)
        _petsc_pchypre_hook(slesp->name, slesp, is_symm, pc);
#else
        cs_base_warn(__FILE__, __LINE__);
        cs_log_printf(CS_LOG_DEFAULT,
                      "%s: Eq. %s: Switch to GAMG since BoomerAMG is not"
                      " available.\n",
                      __func__, slesp->name);
        _petsc_pcgamg_hook(slesp->name, slesp, is_symm, pc);
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

  /* Apply modifications to the PC structure given with command lines.
   * This setting stands for a first setting and may be overwritten with
   * parameters stored in the structure cs_param_sles_t
   * To get the last word use cs_user_sles_petsc_hook()  */

  PCSetFromOptions(pc);
  PCSetUp(pc);
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
  /* No choice otherwise PETSc yields an error */

  slesp->resnorm_type = CS_PARAM_RESNORM_NORM2_RHS;
  KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED);

  /* 2) Set the krylov solver */

  switch (slesp->solver) {

  case CS_PARAM_ITSOL_NONE:
    KSPSetType(ksp, KSPPREONLY);
    break;

  case CS_PARAM_ITSOL_BICG:      /* Improved Bi-CG stab */
    KSPSetType(ksp, KSPIBCGS);
    break;

  case CS_PARAM_ITSOL_BICGSTAB2: /* Preconditioned BiCGstab2 */
    KSPSetType(ksp, KSPBCGSL);
    break;

  case CS_PARAM_ITSOL_CG:        /* Preconditioned Conjugate Gradient */
    if (slesp->precond == CS_PARAM_PRECOND_AMG)
      KSPSetType(ksp, KSPFCG);
    else
      KSPSetType(ksp, KSPCG);
    break;

  case CS_PARAM_ITSOL_FCG:       /* Flexible Conjuguate Gradient */
    KSPSetType(ksp, KSPFCG);
    break;

  case CS_PARAM_ITSOL_FGMRES:    /* Preconditioned flexible GMRES */
    KSPSetType(ksp, KSPFGMRES);
    break;

  case CS_PARAM_ITSOL_GCR:       /* Generalized Conjuguate Residual */
    KSPSetType(ksp, KSPGCR);
    break;

  case CS_PARAM_ITSOL_GMRES:     /* Preconditioned GMRES */
    KSPSetType(ksp, KSPLGMRES);
    break;

  case CS_PARAM_ITSOL_MINRES:    /* Minimal residual */
    KSPSetType(ksp, KSPMINRES);
    break;

  case CS_PARAM_ITSOL_MUMPS:     /* Direct solver (factorization) */
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
    _petsc_cmd(true, slesp->name, "ksp_gmres_modifiedgramschmidt", "1");
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
  case CS_PARAM_ITSOL_FGMRES: /* Flexible GMRES */
    KSPGMRESSetRestart(ksp, slesp->restart);
    break;

  case CS_PARAM_ITSOL_GCR: /* Preconditioned GCR */
    KSPGCRSetRestart(ksp, slesp->restart);
    break;

#if defined(PETSC_HAVE_MUMPS)
  case CS_PARAM_ITSOL_MUMPS:
    {
      cs_param_sles_mumps_t  *mumpsp = slesp->context_param;
      assert(mumpsp != NULL);

      PC  pc;
      KSPGetPC(ksp, &pc);

      if (mumpsp->facto_type == CS_PARAM_SLES_FACTO_LU) {

        PCSetType(pc, PCLU);
        PCFactorSetMatSolverType(pc, MATSOLVERMUMPS);

      }
      else {

        assert(mumpsp->facto_type == CS_PARAM_SLES_FACTO_LDLT_SPD ||
               mumpsp->facto_type == CS_PARAM_SLES_FACTO_LDLT_SYM);

        if (mumpsp->facto_type == CS_PARAM_SLES_FACTO_LDLT_SPD) {

          /* Retrieve the matrices related to this KSP */

          Mat a, pa;
          KSPGetOperators(ksp, &a, &pa);
          MatSetOption(a, MAT_SPD, PETSC_TRUE); /* set MUMPS id%SYM=1 */

        }

        PCSetType(pc, PCCHOLESKY);
        PCFactorSetMatSolverType(pc, MATSOLVERMUMPS);
        PCFactorSetUpMatSolverType(pc); /* call MatGetFactor() to create F */

      } /* L.D.Lt factorization */
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
                   slesp->cvg_param.rtol,  /* relative convergence tolerance */
                   slesp->cvg_param.atol,  /* absolute convergence tolerance */
                   slesp->cvg_param.dtol,  /* divergence tolerance */
                   slesp->cvg_param.n_max_iter);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set PETSc solver and preconditioner
 *
 * \param[in, out] context    pointer to optional (untyped) value or structure
 * \param[in, out] ksp_struct pointer to PETSc KSP context
 */
/*----------------------------------------------------------------------------*/

static void
_petsc_setup_hook(void    *context,
                  void    *ksp_struct)
{
  cs_param_sles_t  *slesp = (cs_param_sles_t  *)context;
  KSP  ksp = ksp_struct;

  cs_fp_exception_disable_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */

  int len = strlen(slesp->name) + 1;
  char  *prefix = NULL;
  BFT_MALLOC(prefix, len + 1, char);
  sprintf(prefix, "%s_", slesp->name);
  prefix[len] = '\0';
  KSPSetOptionsPrefix(ksp, prefix);
  BFT_FREE(prefix);

  /* 1) Set the solver */

  _petsc_set_krylov_solver(slesp, ksp);

  /* 2) Set the preconditioner */

  _petsc_set_pc_type(slesp, ksp);

  /* 3) User function for additional settings */

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
 * \param[in]      slesp       pointer to the SLES parameter settings
 * \param[in, out] ksp         pointer to PETSc KSP structure (solver)
 */
/*----------------------------------------------------------------------------*/

static void
_petsc_common_block_hook(const cs_param_sles_t    *slesp,
                         KSP                       ksp)
{
  PC  pc;
  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCFIELDSPLIT);

  switch (slesp->pcd_block_type) {
  case CS_PARAM_PRECOND_BLOCK_UPPER_TRIANGULAR:
  case CS_PARAM_PRECOND_BLOCK_LOWER_TRIANGULAR:
  case CS_PARAM_PRECOND_BLOCK_FULL_UPPER_TRIANGULAR:
  case CS_PARAM_PRECOND_BLOCK_FULL_LOWER_TRIANGULAR:
    PCFieldSplitSetType(pc, PC_COMPOSITE_MULTIPLICATIVE);
    break;

  case CS_PARAM_PRECOND_BLOCK_SYM_GAUSS_SEIDEL:
  case CS_PARAM_PRECOND_BLOCK_FULL_SYM_GAUSS_SEIDEL:
    PCFieldSplitSetType(pc, PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE);
    break;

  case CS_PARAM_PRECOND_BLOCK_DIAG:
  case CS_PARAM_PRECOND_BLOCK_FULL_DIAG:
  default:
    PCFieldSplitSetType(pc, PC_COMPOSITE_ADDITIVE);
    break;
  }

  /* Apply modifications to the KSP structure */

  PCFieldSplitSetBlockSize(pc, 3);

  PetscInt  id = 0;
  PCFieldSplitSetFields(pc, "x", 1, &id, &id);
  id = 1;
  PCFieldSplitSetFields(pc, "y", 1, &id, &id);
  id = 2;
  PCFieldSplitSetFields(pc, "z", 1, &id, &id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer: setup hook for setting PETSc solver and
 *         preconditioner.
 *         Case of multiplicative AMG block preconditioner for a CG with GAMG
 *         as AMG type
 *
 * \param[in, out] context    pointer to optional (untyped) value or structure
 * \param[in, out] ksp_struct pointer to PETSc KSP context
 */
/*----------------------------------------------------------------------------*/

static void
_petsc_amg_block_gamg_hook(void     *context,
                           void     *ksp_struct)
{
  cs_param_sles_t  *slesp = (cs_param_sles_t *)context;
  KSP  ksp = ksp_struct;

  cs_fp_exception_disable_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */

  /* prefix will be extended with the fieldsplit */

  int len = strlen(slesp->name) + 1;
  int _len = len + strlen("_fieldsplit_x_") + 1;
  char  *prefix = NULL;
  BFT_MALLOC(prefix, _len + 1, char);
  sprintf(prefix, "%s_", slesp->name);
  prefix[len] = '\0';
  KSPSetOptionsPrefix(ksp, prefix);

  /* Set the solver */

  _petsc_set_krylov_solver(slesp, ksp);

  /* Common settings to block preconditionner */

  _petsc_common_block_hook(slesp, ksp);

  PC  pc;
  KSPGetPC(ksp, &pc);
  PCSetUp(pc);

  PC  _pc;
  PetscInt  n_split;
  KSP  *xyz_subksp;
  const char xyz[3] = "xyz";
  bool  is_symm = _system_should_be_sym(slesp->solver);

  PCFieldSplitGetSubKSP(pc, &n_split, &xyz_subksp);
  assert(n_split == 3);

  for (PetscInt id = 0; id < n_split; id++) {

    sprintf(prefix, "%s_fieldsplit_%c", slesp->name, xyz[id]);

    _petsc_cmd(true, prefix, "ksp_type","preonly");

    /* Predefined settings when using AMG as a preconditioner */

    KSP  _ksp = xyz_subksp[id];
    KSPGetPC(_ksp, &_pc);

    _petsc_pcgamg_hook(prefix, slesp, is_symm, _pc);

    PCSetFromOptions(_pc);
    KSPSetFromOptions(_ksp);

  } /* Loop on block settings */

  BFT_FREE(prefix);
  PetscFree(xyz_subksp);

  /* User function for additional settings */

  cs_user_sles_petsc_hook(context, ksp);

  PCSetFromOptions(pc);
  KSPSetFromOptions(ksp);

  /* Dump the setup related to PETSc in a specific file */

  if (!slesp->setup_done) {
    KSPSetUp(ksp);
    cs_sles_petsc_log_setup(ksp);
    slesp->setup_done = true;
  }

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
 * \param[in, out] context     pointer to optional (untyped) value or structure
 * \param[in, out] ksp_struct  pointer to PETSc KSP context
 */
/*----------------------------------------------------------------------------*/

static void
_petsc_amg_block_boomer_hook(void     *context,
                             void     *ksp_struct)
{
  cs_param_sles_t  *slesp = (cs_param_sles_t *)context;
  KSP  ksp = ksp_struct;

  cs_fp_exception_disable_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */

  /* prefix will be extended with the fieldsplit */

  int len = strlen(slesp->name) + 1;
  int _len = len + strlen("_fieldsplit_x_") + 1;
  char  *prefix = NULL;
  BFT_MALLOC(prefix, _len + 1, char);
  sprintf(prefix, "%s_", slesp->name);
  prefix[len] = '\0';
  KSPSetOptionsPrefix(ksp, prefix);

  /* Set the solver */

  _petsc_set_krylov_solver(slesp, ksp);

  /* Common settings to block preconditionner */

  _petsc_common_block_hook(slesp, ksp);

  /* Predefined settings when using AMG as a preconditioner */

  PC  pc;
  KSPGetPC(ksp, &pc);
  PCSetUp(pc);

  PC  _pc;
  PetscInt  n_split;
  KSP  *xyz_subksp;
  const char xyz[3] = "xyz";
  bool  is_symm = _system_should_be_sym(slesp->solver);

  PCFieldSplitGetSubKSP(pc, &n_split, &xyz_subksp);
  assert(n_split == 3);

  for (PetscInt id = 0; id < n_split; id++) {

    sprintf(prefix, "%s_fieldsplit_%c", slesp->name, xyz[id]);

    _petsc_cmd(true, prefix, "ksp_type", "preonly");

    /* Predefined settings when using AMG as a preconditioner */

    KSP  _ksp = xyz_subksp[id];
    KSPGetPC(_ksp, &_pc);

    _petsc_pchypre_hook(prefix, slesp, is_symm, _pc);

    PCSetFromOptions(_pc);
    KSPSetFromOptions(_ksp);

  } /* Loop on block settings */

  BFT_FREE(prefix);
  PetscFree(xyz_subksp);

  /* User function for additional settings */

  cs_user_sles_petsc_hook(context, ksp);

  PCSetFromOptions(pc);
  KSPSetFromOptions(ksp);

  /* Dump the setup related to PETSc in a specific file */

  if (!slesp->setup_done) {
    KSPSetUp(ksp);
    cs_sles_petsc_log_setup(ksp);
    slesp->setup_done = true;
  }

  cs_fp_exception_restore_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */
}
#endif

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer: setup hook for setting PETSc solver and
 *         preconditioner.
 *         Case of block preconditioner
 *
 * \param[in, out] context    pointer to optional (untyped) value or structure
 * \param[in, out] ksp_struct pointer to PETSc KSP context
 */
/*----------------------------------------------------------------------------*/

static void
_petsc_block_hook(void     *context,
                  void     *ksp_struct)
{
  cs_param_sles_t  *slesp = (cs_param_sles_t *)context;
  KSP  ksp = ksp_struct;

  cs_fp_exception_disable_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */

  /* prefix will be extended with the fieldsplit */

  int len = strlen(slesp->name) + 1;
  int _len = len + strlen("_fieldsplit_x_") + 1;
  char  *prefix = NULL;
  BFT_MALLOC(prefix, _len + 1, char);
  sprintf(prefix, "%s_", slesp->name);
  prefix[len] = '\0';
  KSPSetOptionsPrefix(ksp, prefix);

  /* Set the solver (tolerance and max_it too) */

  _petsc_set_krylov_solver(slesp, ksp);

  /* Common settings to block preconditionner */

  _petsc_common_block_hook(slesp, ksp);

  PC  pc;
  KSPGetPC(ksp, &pc);
  PCSetUp(pc);

  PC  _pc;
  KSP  *xyz_subksp;
  PetscInt  n_split;
  const char  xyz[3] = "xyz";

  PCFieldSplitGetSubKSP(pc, &n_split, &xyz_subksp);
  assert(n_split == 3);

  for (PetscInt id = 0; id < n_split; id++) {

    sprintf(prefix, "%s_fieldsplit_%c", slesp->name, xyz[id]);

    switch (slesp->precond) {

    case CS_PARAM_PRECOND_NONE:
      _petsc_cmd(true, prefix, "ksp_type", "richardson");
      break;

    case CS_PARAM_PRECOND_DIAG:
      _petsc_cmd(true, prefix, "ksp_type", "richardson");
      _petsc_cmd(true, prefix, "pc_type", "jacobi");
      break;

    case CS_PARAM_PRECOND_ILU0:
    case CS_PARAM_PRECOND_BJACOB_ILU0:
      if (slesp->solver_class == CS_PARAM_SLES_CLASS_HYPRE) {
#if defined(PETSC_HAVE_HYPRE)
        _petsc_cmd(true, prefix, "ksp_type", "preonly");
        _petsc_cmd(true, prefix, "pc_type", "hypre");
        _petsc_cmd(true, prefix, "pc_hypre_type", "euclid");
        _petsc_cmd(true, prefix, "pc_hypre_euclid_level", "0");
#else
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid option: HYPRE is not installed.", __func__);
#endif
      }
      else {

          _petsc_cmd(true, prefix, "ksp_type", "richardson");
          _petsc_bilu0_hook(prefix);

      }
      break;

    case CS_PARAM_PRECOND_ICC0:
      _petsc_cmd(true, prefix, "ksp_type", "richardson");
      _petsc_bicc0_hook(prefix);
      break;

    case CS_PARAM_PRECOND_LU:
      _petsc_cmd(true, prefix, "ksp_type", "preonly");
#if defined(PETSC_HAVE_MUMPS)
      _petsc_cmd(true, prefix, "pc_type", "lu");
      _petsc_cmd(true, prefix, "pc_factor_mat_solver_type", "mumps");
#else
      if (cs_glob_n_ranks == 1)
        _petsc_cmd(true, prefix, "pc_type", "lu");
      else { /* Switch to a block version (sequential in each block) */
        _petsc_cmd(true, prefix, "pc_type", "bjacobi");
        _petsc_cmd(true, prefix, "pc_jacobi_blocks", "1");
        _petsc_cmd(true, prefix, "sub_ksp_type", "preonly");
        _petsc_cmd(true, prefix, "sub_pc_type", "lu");
      }
#endif
      break;

    case CS_PARAM_PRECOND_SSOR:
    case CS_PARAM_PRECOND_BJACOB_SGS:
      _petsc_cmd(true, prefix, "ksp_type", "richardson");
      _petsc_bssor_hook(prefix);
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " %s: Eq. %s: Invalid preconditioner.",__func__, slesp->name);
      break;

    } /* Switch on preconditioning type */

    KSP  _ksp = xyz_subksp[id];
    KSPGetPC(_ksp, &_pc);
    PCSetFromOptions(_pc);
    KSPSetUp(_ksp);

  } /* Loop on block settings */

  BFT_FREE(prefix);
  PetscFree(xyz_subksp);

  /* User function for additional settings */

  cs_user_sles_petsc_hook(context, ksp);

  PCSetFromOptions(pc);
  KSPSetFromOptions(ksp);

  /* Dump the setup related to PETSc in a specific file */

  if (!slesp->setup_done) {
    KSPSetUp(ksp);
    cs_sles_petsc_log_setup(ksp);
    slesp->setup_done = true;
  }

  cs_fp_exception_restore_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */
}
#endif /* defined(HAVE_PETSC) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check if the settings are consistent. Can apply minor modifications.
 *
 * \param[in, out]  slesp       pointer to a \ref cs_param_sles_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_check_settings(cs_param_sles_t     *slesp)
{
  if (cs_param_sles_is_mumps_set(slesp->solver)) { /* Checks related to MUMPS */

    cs_param_sles_class_t  ret_class =
      cs_param_sles_check_class(CS_PARAM_SLES_CLASS_MUMPS);
    if (ret_class == CS_PARAM_SLES_N_CLASSES)
      bft_error(__FILE__, __LINE__, 0,
                " %s: Error detected while setting the SLES \"%s\"\n"
                " MUMPS is not available with your installation.\n"
                " Please check your installation settings.\n",
                __func__, slesp->name);
    else
      slesp->solver_class = ret_class;

  }
  else {
    if (slesp->solver_class == CS_PARAM_SLES_CLASS_MUMPS)
      bft_error(__FILE__, __LINE__, 0,
                " %s: Error detected while setting the SLES \"%s\"\n"
                " MUMPS class is not consistent with your settings.\n"
                " Please check your installation settings.\n",
                __func__, slesp->name);
  }

  /* Checks related to GCR/GMRES algorithms */

  if (slesp->solver == CS_PARAM_ITSOL_GMRES ||
      slesp->solver == CS_PARAM_ITSOL_FGMRES ||
      slesp->solver == CS_PARAM_ITSOL_GCR)
    if (slesp->restart < 2)
      bft_error(__FILE__, __LINE__, 0,
                " %s: Error detected while setting the SLES \"%s\"\n"
                " The restart interval (=%d) is not big enough.\n"
                " Please check your installation settings.\n",
                __func__, slesp->name, slesp->restart);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the value of the polynomial degree to consider according to
 *        the settings. Only for in-house solvers.
 *
 * \param[in]  slesp      pointer to a \ref cs_param_sles_t structure
 *
 * \return the value of the polynomial degree to consider
 */
/*----------------------------------------------------------------------------*/

static int
_get_poly_degree(const cs_param_sles_t     *slesp)
{
  switch (slesp->precond) { /* Switch on the preconditioner type */

  case CS_PARAM_PRECOND_DIAG:
    return 0;

  case CS_PARAM_PRECOND_POLY1:
    return 1;

  case CS_PARAM_PRECOND_POLY2:
    return 2;

  default:
    return -1;

  }
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

static void
_set_saturne_sles(bool                 use_field_id,
                  cs_param_sles_t     *slesp)
{
  assert(slesp != NULL);  /* Sanity checks */

  const char  *sles_name = use_field_id ? NULL : slesp->name;
  assert(slesp->field_id > -1 || sles_name != NULL);

  /* Retrieve the sles structure for this equation */

  cs_sles_t  *sles = cs_sles_find(slesp->field_id, sles_name);

  if (sles != NULL) {
    if (slesp->field_id > -1) {
      /* Solver settings already forced */
      return;
    }
  }
  else
    sles = cs_sles_find_or_add(slesp->field_id, sles_name);

  int  poly_degree = _get_poly_degree(slesp);

  /* Retrieve associated context structures */

  cs_sles_it_t  *itsol = cs_sles_get_context(sles);
  cs_sles_pc_t  *pc = cs_sles_it_get_pc(itsol);
  cs_multigrid_t  *mg = NULL;

  /* 1- Define the iterative solver
   *    =========================== */

  if (itsol == NULL) { /* Not already defined. Add a new one. */

    switch (slesp->solver) {

    case CS_PARAM_ITSOL_AMG:
      {
        switch (slesp->amg_type) {

        case CS_PARAM_AMG_HOUSE_V:
          mg = cs_multigrid_define(slesp->field_id, sles_name,
                                   CS_MULTIGRID_V_CYCLE);

          /* Advanced setup (default is specified inside the brackets)
           * for AMG as solver */

          cs_multigrid_set_solver_options
            (mg,
             CS_SLES_JACOBI,              /* descent smoother (CS_SLES_PCG) */
             CS_SLES_JACOBI,              /* ascent smoother (CS_SLES_PCG) */
             CS_SLES_PCG,                 /* coarse solver (CS_SLES_PCG) */
             slesp->cvg_param.n_max_iter, /* n max cycles (100) */
             5,       /* n max iter for descent (10) */
             5,       /* n max iter for ascent (10) */
             1000,    /* n max iter coarse solver (10000) */
             0,       /* polynomial precond. degree descent (0) */
             0,       /* polynomial precond. degree ascent (0) */
             -1,      /* polynomial precond. degree coarse (0) */
             1.0,     /* precision multiplier descent (< 0 forces max iters) */
             1.0,     /* precision multiplier ascent (< 0 forces max iters) */
             1);      /* requested precision multiplier coarse (default 1) */
          break;

        case CS_PARAM_AMG_HOUSE_K:
          mg = cs_multigrid_define(slesp->field_id, sles_name,
                                   CS_MULTIGRID_K_CYCLE);

          cs_multigrid_set_solver_options
            (mg,
             CS_SLES_P_SYM_GAUSS_SEIDEL,    /* descent smoother */
             CS_SLES_P_SYM_GAUSS_SEIDEL,    /* ascent smoother */
             CS_SLES_PCG,                   /* coarse smoother */
             slesp->cvg_param.n_max_iter,   /* n_max_cycles */
             1,                             /* n_max_iter_descent, */
             1,                             /* n_max_iter_ascent */
             100,                           /* n_max_iter_coarse */
             0,                             /* poly_degree_descent */
             0,                             /* poly_degree_ascent */
             0,                             /* poly_degree_coarse */
             -1.0,                          /* precision_mult_descent */
             -1.0,                          /* precision_mult_ascent */
             1);                            /* precision_mult_coarse */
          break;

        default:
          bft_error(__FILE__, __LINE__, 0,
                    " %s; System: %s -- Invalid AMG type with code_saturne"
                    " solvers.", __func__, slesp->name);
          break;

        } /* End of switch on the AMG type */

      }
      break; /* AMG as solver */

    case CS_PARAM_ITSOL_BICG:
      itsol = cs_sles_it_define(slesp->field_id, sles_name,
                                CS_SLES_BICGSTAB,
                                poly_degree,
                                slesp->cvg_param.n_max_iter);
      break;

    case CS_PARAM_ITSOL_BICGSTAB2:
      itsol = cs_sles_it_define(slesp->field_id, sles_name,
                                CS_SLES_BICGSTAB2,
                                poly_degree,
                                slesp->cvg_param.n_max_iter);
      break;

    case CS_PARAM_ITSOL_CG:
      if (slesp->flexible) {
        slesp->solver = CS_PARAM_ITSOL_FCG;
        itsol = cs_sles_it_define(slesp->field_id, sles_name,
                                  CS_SLES_IPCG,
                                  poly_degree,
                                  slesp->cvg_param.n_max_iter);
      }
      else
        itsol = cs_sles_it_define(slesp->field_id, sles_name,
                                  CS_SLES_PCG,
                                  poly_degree,
                                  slesp->cvg_param.n_max_iter);
      break;

    case CS_PARAM_ITSOL_CR3:
      itsol = cs_sles_it_define(slesp->field_id, sles_name,
                                CS_SLES_PCR3,
                                poly_degree,
                                slesp->cvg_param.n_max_iter);
      break;

    case CS_PARAM_ITSOL_FCG:
      itsol = cs_sles_it_define(slesp->field_id, sles_name,
                                CS_SLES_FCG,
                                poly_degree,
                                slesp->cvg_param.n_max_iter);
      break;

    case CS_PARAM_ITSOL_GKB_CG:
      itsol = cs_sles_it_define(slesp->field_id, sles_name,
                                CS_SLES_IPCG,
                                poly_degree,
                                slesp->cvg_param.n_max_iter);
      break;

    case CS_PARAM_ITSOL_GAUSS_SEIDEL:
      itsol = cs_sles_it_define(slesp->field_id, sles_name,
                                CS_SLES_P_GAUSS_SEIDEL,
                                -1, /* Not useful to apply a preconditioner */
                                slesp->cvg_param.n_max_iter);
      break;

    case CS_PARAM_ITSOL_FGMRES:  /* Not available --> close to GCR */
      cs_base_warn(__FILE__, __LINE__);
      bft_printf(" Switch to the GCR implementation of code_saturne\n");
      /* No break (wanted behavior) */
    case CS_PARAM_ITSOL_GKB_GMRES:
    case CS_PARAM_ITSOL_GCR:
      itsol = cs_sles_it_define(slesp->field_id, sles_name,
                                CS_SLES_GCR,
                                poly_degree,
                                slesp->cvg_param.n_max_iter);
      break;

    case CS_PARAM_ITSOL_GMRES:
      if (slesp->flexible) {
        slesp->solver = CS_PARAM_ITSOL_GCR;
        itsol = cs_sles_it_define(slesp->field_id, sles_name,
                                  CS_SLES_GCR,
                                  poly_degree,
                                  slesp->cvg_param.n_max_iter);
      }
      else
        itsol = cs_sles_it_define(slesp->field_id, sles_name,
                                  CS_SLES_GMRES,
                                  poly_degree,
                                  slesp->cvg_param.n_max_iter);
      break;

    case CS_PARAM_ITSOL_JACOBI:
      itsol = cs_sles_it_define(slesp->field_id, sles_name,
                                CS_SLES_JACOBI,
                                -1, /* Not useful to apply a preconditioner */
                                slesp->cvg_param.n_max_iter);
      break;

    case CS_PARAM_ITSOL_SYM_GAUSS_SEIDEL:
      itsol = cs_sles_it_define(slesp->field_id, sles_name,
                                CS_SLES_P_SYM_GAUSS_SEIDEL,
                                -1, /* Not useful to apply a preconditioner */
                                slesp->cvg_param.n_max_iter);
      break;

    case CS_PARAM_ITSOL_USER_DEFINED:
      itsol = cs_sles_it_define(slesp->field_id, sles_name,
                                CS_SLES_USER_DEFINED,
                                poly_degree,
                                slesp->cvg_param.n_max_iter);
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid iterative solver for solving equation %s.\n"
                " Please modify your settings.", __func__, slesp->name);
      break;

    } /* End of switch */

  }
  else {

    /* The itsol structure has already been defined. Check that this is
       consistent with the sles parameters */

    cs_sles_it_type_t  itsol_type = cs_sles_it_get_type(itsol);
    int  ret_code = -1;
    switch (itsol_type) {

    case CS_SLES_PCG:
      if (slesp->solver != CS_PARAM_ITSOL_CG)
        ret_code = 0;
      break;

    case CS_SLES_FCG:
    case CS_SLES_IPCG:
      if (slesp->solver != CS_PARAM_ITSOL_FCG &&
          slesp->solver != CS_PARAM_ITSOL_GKB_CG)
        ret_code = 1;
      break;

    case CS_SLES_JACOBI:
      if (slesp->solver != CS_PARAM_ITSOL_JACOBI)
        ret_code = 2;
      break;

    case CS_SLES_BICGSTAB:
      if (slesp->solver != CS_PARAM_ITSOL_BICG)
        ret_code = 3;
      break;

    case CS_SLES_BICGSTAB2:
      if (slesp->solver != CS_PARAM_ITSOL_BICGSTAB2)
        ret_code = 4;
      break;

    case CS_SLES_GCR:
      if (slesp->solver != CS_PARAM_ITSOL_GCR &&
          slesp->solver != CS_PARAM_ITSOL_GKB_GMRES)
        ret_code = 5;
      break;

    case CS_SLES_GMRES:
      if (slesp->solver != CS_PARAM_ITSOL_GMRES)
        ret_code = 6;
      break;

    case CS_SLES_P_GAUSS_SEIDEL:
      if (slesp->solver != CS_PARAM_ITSOL_GAUSS_SEIDEL)
        ret_code = 7;
      break;

    case CS_SLES_P_SYM_GAUSS_SEIDEL:
      if (slesp->solver != CS_PARAM_ITSOL_SYM_GAUSS_SEIDEL)
        ret_code = 8;
      break;

    case CS_SLES_PCR3:
      if (slesp->solver != CS_PARAM_ITSOL_CR3)
        ret_code = 9;
      break;

    case CS_SLES_USER_DEFINED:
      if (slesp->solver != CS_PARAM_ITSOL_USER_DEFINED)
        ret_code = 10;
      break;

    default:
      ret_code = 11;
      break;

    } /* End of switch on itsol_type */

    if (ret_code > -1)
      bft_error(__FILE__, __LINE__, 0,
                "%s: Invalid solver w.r.t. settings (code: %d)\n",
                __func__, ret_code);

  } /* Iterative solver already defined */

  if (slesp->flexible) { /* Additional checks */

    switch (cs_sles_it_get_type(itsol)) {

    case CS_SLES_PCG:
    case CS_SLES_BICGSTAB:
    case CS_SLES_BICGSTAB2:
    case CS_SLES_GMRES:
    case CS_SLES_PCR3:
      cs_base_warn(__FILE__, __LINE__);
      bft_printf("--> A flexible Krylov method should be used.\n");
      break;

    default:
      break;
    }

  }

  /* 2- Define the preconditioner
   *    ========================= */

  if (pc == NULL) { /* Add the preconditioner context */

    assert(itsol != NULL);

    switch (slesp->precond) {

    case CS_PARAM_PRECOND_AMG:
      /* -------------------- */
      if (slesp->amg_type == CS_PARAM_AMG_HOUSE_V) {

        pc = cs_multigrid_pc_create(CS_MULTIGRID_V_CYCLE);
        mg = cs_sles_pc_get_context(pc);
        cs_sles_it_transfer_pc(itsol, &pc);

      }
      else if (slesp->amg_type == CS_PARAM_AMG_HOUSE_K) {

        pc = cs_multigrid_pc_create(CS_MULTIGRID_K_CYCLE);
        mg = cs_sles_pc_get_context(pc);
        cs_sles_it_transfer_pc(itsol, &pc);

        /* Change the default settings when used as preconditioner */

        cs_multigrid_set_solver_options
          (mg,
           CS_SLES_PCG,   /* descent smoother */
           CS_SLES_PCG,   /* ascent smoother */
           CS_SLES_PCG,   /* coarse solver */
           1,             /* n_max_cycles */
           2,             /* n_max_iter_descent, */
           2,             /* n_max_iter_ascent */
           500,           /* n_max_iter_coarse */
           0,             /* poly_degree_descent */
           0,             /* poly_degree_ascent */
           0,             /* poly_degree_coarse */
           -1.0,          /* precision_mult_descent */
           -1.0,          /* precision_mult_ascent */
           1.0);          /* precision_mult_coarse */

        cs_multigrid_set_coarsening_options(mg,
                                            8,    /* aggregation_limit */
                                            CS_GRID_COARSENING_SPD_PW,
                                            10,   /* n_max_levels */
                                            150,  /* min_g_cells */
                                            0.,   /* P0P1 relaxation */
                                            0);   /* postprocess */

      }
      else
        bft_error(__FILE__, __LINE__, 0, " %s: System: %s;"
                  " Invalid AMG type with code_saturne solvers.",
                  __func__, slesp->name);
      break;

    case CS_PARAM_PRECOND_MUMPS:
#if defined(HAVE_MUMPS)
      if (slesp->context_param == NULL)
        cs_param_sles_mumps(slesp,
                            true, /* single by default in case if precond. */
                            CS_PARAM_SLES_FACTO_LU);

      pc = cs_sles_mumps_pc_create(slesp);
#else
      bft_error(__FILE__, __LINE__, 0, "%s: MUMPS not available in this build.",
                __func__);
#endif
      cs_sles_it_transfer_pc(itsol, &pc);
      break;

    default: /* Nothing else to do */
      break;

    } /* Switch on the preconditioner */

  } /* Preconditioner is not defined */

  /* In case of high verbosity, additional output are generated */

  if (slesp->verbosity > 3) /* true=use_iteration instead of wall clock time */
    cs_sles_it_set_plot_options(itsol, slesp->name, true);
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

static void
_set_mumps_sles(bool                 use_field_id,
                cs_param_sles_t     *slesp)
{
  assert(slesp != NULL);  /* Sanity checks */

  const  char  *sles_name = use_field_id ? NULL : slesp->name;
  assert(slesp->field_id > -1 || sles_name != NULL);

  if (slesp->context_param == NULL) /* Define the context by default */
    cs_param_sles_mumps(slesp,
                        false,  /* is single-precision ? */
                        CS_PARAM_SLES_FACTO_LU);

#if defined(HAVE_MUMPS)
  cs_sles_mumps_define(slesp->field_id,
                       sles_name,
                       slesp,
                       cs_user_sles_mumps_hook,
                       NULL);
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

static void
_set_petsc_hypre_sles(bool                 use_field_id,
                      cs_param_sles_t     *slesp)
{
  assert(slesp != NULL);  /* Sanity checks */

  const  char  *sles_name = use_field_id ? NULL : slesp->name;
  assert(slesp->field_id > -1 || sles_name != NULL);

#if defined(HAVE_PETSC)
  cs_sles_petsc_init();

  if (slesp->pcd_block_type != CS_PARAM_PRECOND_BLOCK_NONE) {

    if (slesp->precond == CS_PARAM_PRECOND_AMG) {

      if (slesp->amg_type == CS_PARAM_AMG_PETSC_GAMG_V ||
          slesp->amg_type == CS_PARAM_AMG_PETSC_GAMG_W)
        cs_sles_petsc_define(slesp->field_id,
                             sles_name,
                             MATMPIAIJ,
                             _petsc_amg_block_gamg_hook,
                             (void *)slesp);

      else if (slesp->amg_type == CS_PARAM_AMG_HYPRE_BOOMER_V ||
               slesp->amg_type == CS_PARAM_AMG_HYPRE_BOOMER_W) {
#if defined(PETSC_HAVE_HYPRE)
        cs_sles_petsc_define(slesp->field_id,
                             sles_name,
                             MATMPIAIJ,
                             _petsc_amg_block_boomer_hook,
                             (void *)slesp);
#else
        cs_base_warn(__FILE__, __LINE__);
        cs_log_printf(CS_LOG_DEFAULT,
                      " %s: System: %s.\n"
                      " Boomer is not available. Switch to GAMG solver.",
                      __func__, slesp->name);
        cs_sles_petsc_define(slesp->field_id,
                             sles_name,
                             MATMPIAIJ,
                             _petsc_amg_block_gamg_hook,
                             (void *)slesp);
#endif
      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  " %s: System: %s\n"
                  " No AMG solver available for a block-AMG.",
                  __func__, slesp->name);
    }
    else {

      cs_sles_petsc_define(slesp->field_id,
                           sles_name,
                           MATMPIAIJ,
                           _petsc_block_hook,
                           (void *)slesp);

    }

  }
  else { /* No block preconditioner */

    cs_sles_petsc_define(slesp->field_id,
                         sles_name,
                         MATMPIAIJ,
                         _petsc_setup_hook,
                         (void *)slesp);

  }
#else
  bft_error(__FILE__, __LINE__, 0,
            _(" %s: PETSC algorithms used to solve %s are not linked.\n"
              " Please install code_saturne with PETSc."),
            __func__, slesp->name);
#endif /* HAVE_PETSC */
}

#if defined(HAVE_HYPRE)
/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the solver when HYPRE is used
 *         Check HYPRE documentation for available options:
 *         https://hypre.readthedocs.io/en/latest/index.html
 *
 * \param[in]      slesp   pointer to a set of SLES parameters
 * \param[in, out] hs      pointer to the HYPRE solver structure
 *
 * \return a HYPRE_Solver structure (pointer) for the preconditioner
 */
/*----------------------------------------------------------------------------*/

static HYPRE_Solver
_set_hypre_solver(cs_param_sles_t    *slesp,
                  HYPRE_Solver        hs)
{
  HYPRE_Solver  pc = NULL;

  switch (slesp->solver) {

  case CS_PARAM_ITSOL_AMG: /* BoomerAMG as solver. Nothing to do at this
                              stage. This is done in the calling function. */
    pc = hs;
    break;

  case CS_PARAM_ITSOL_BICG:
  case CS_PARAM_ITSOL_BICGSTAB2:
    HYPRE_BiCGSTABSetTol(hs, (HYPRE_Real)slesp->cvg_param.rtol);
    HYPRE_BiCGSTABSetMaxIter(hs, (HYPRE_Int)slesp->cvg_param.n_max_iter);
    HYPRE_BiCGSTABGetPrecond(hs, &pc);
    break;

  case CS_PARAM_ITSOL_CG:
  case CS_PARAM_ITSOL_FCG:
    HYPRE_PCGSetTol(hs, (HYPRE_Real)slesp->cvg_param.rtol);
    HYPRE_PCGSetMaxIter(hs, (HYPRE_Int)slesp->cvg_param.n_max_iter);
    HYPRE_PCGGetPrecond(hs, &pc);
    break;

  case CS_PARAM_ITSOL_FGMRES:
  case CS_PARAM_ITSOL_GCR:
    HYPRE_FlexGMRESSetTol(hs, (HYPRE_Real)slesp->cvg_param.rtol);
    HYPRE_FlexGMRESSetMaxIter(hs, (HYPRE_Int)slesp->cvg_param.n_max_iter);
    HYPRE_FlexGMRESSetKDim(hs, (HYPRE_Int)slesp->restart);
    HYPRE_FlexGMRESGetPrecond(hs, &pc);
    break;

  case CS_PARAM_ITSOL_GMRES:
    HYPRE_GMRESSetTol(hs, (HYPRE_Real)slesp->cvg_param.rtol);
    HYPRE_GMRESSetMaxIter(hs, (HYPRE_Int)slesp->cvg_param.n_max_iter);
    HYPRE_GMRESSetKDim(hs, (HYPRE_Int)slesp->restart);
    HYPRE_GMRESGetPrecond(hs, &pc);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid type of solver for eq. \"%s\"\n",
              __func__, slesp->name);
    break;
  }

  return pc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup hook function for a Hypre KSP solver with preconditioner
 *         This function is called at the end of the setup stage for a KSP
 *         solver.
 *
 *         Note: if the context pointer is non-NULL, it must point to valid
 *         data when the selection function is called so that value or
 *         structure should not be temporary (i.e. local);
 *
 *         Check HYPRE documentation for available options:
 *         https://hypre.readthedocs.io/en/latest/index.html
 *
 * \param[in]      verbosity  verbosity level
 * \param[in, out] context    pointer to optional (untyped) value or structure
 * \param[in, out] solver_p   pointer to a HYPRE solver (cast ont-the-fly)
 */
/*----------------------------------------------------------------------------*/

static void
_hypre_boomeramg_hook(int     verbosity,
                      void   *context,
                      void   *solver_p)
{
  cs_param_sles_t  *slesp = (cs_param_sles_t *)context;
  HYPRE_Solver  hs = solver_p;
  HYPRE_Solver  amg = _set_hypre_solver(slesp, hs);
  bool  amg_as_precond = (slesp->solver == CS_PARAM_ITSOL_AMG) ? false : true;

  /* Set the verbosity level */

  HYPRE_BoomerAMGSetPrintLevel(amg, (verbosity > 3 ? 3 : verbosity));

  /* Set the type of cycle */

  switch (slesp->amg_type) {

  case CS_PARAM_AMG_HYPRE_BOOMER_V:
    HYPRE_BoomerAMGSetCycleType(amg, 1);
    break;

  case CS_PARAM_AMG_HYPRE_BOOMER_W:
    HYPRE_BoomerAMGSetCycleType(amg, 2);
    break;

    /* TODO: Add a full multigrid */

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid type of AMG cycle for eq. \"%s\"\n",
              __func__, slesp->name);
    break;
  }

  /* From HYPRE website. Coarsening type. Default = 10
   *
   * 22 and 8 are also good choices
   *
   *  HMIS (=10): uses one pass Ruge-Stueben on each processor independently,
   *  followed by PMIS using the interior C-points generated as its first
   *  independent set
   *
   *  PMIS (=8) can be also a good choice
   *  a parallel coarsening algorithm using independent sets, generating lower
   *  complexities than CLJP (=0), might also lead to slower convergence
   *
   * 21: CGC coarsening by M. Griebel, B. Metsch and A. Schweitzer
   * 22: CGC-E coarsening by M. Griebel, B. Metsch and A.Schweitzer
   */

  HYPRE_BoomerAMGSetCoarsenType(amg, 10);

  /* From HYPRE documentation.
   * InterpType: recommended: 6 (default) or 7
   *
   *  2: classical modified interpolation for hyperbolic PDEs
   *  4: multipass interpolation
   *  6: extended+i interpolation (also for GPU use)
   *  7: extended+i (if no common C neighbor) interpolation
   * 10: classical block interpolation (with nodal systems only)
   * 11: classical block interpolation (with nodal systems only) + diagonalized
   *     diagonal blocks
   * 13: FF1 interpolation
   * 16: extended interpolation in matrix-matrix form
   * 17: extended+i interpolation in matrix-matrix form
   * 18: extended+e interpolation in matrix-matrix form
   *
   * PMAxElmts: maximal number of elements per row for the interpolation
   *  recommended 4 (default) or 5 --> test with CDO --> 8
   */

  HYPRE_BoomerAMGSetInterpType(amg, 4);
  HYPRE_BoomerAMGSetPMaxElmts(amg, 8);

  /* From the HYPRE documentation.
   * "One important parameter is the strong threshold, which can be set using
   * the function HYPRE_BoomerAMGSetStrongThreshold. The default value is 0.25,
   * which appears to be a good choice for 2- dimensional problems and the low
   * complexity coarsening algorithms. For 3-dimensional problems a better
   * choice appears to be 0.5, when using the default coarsening algorithm.
   * However, the choice of the strength threshold is problem dependent and
   * therefore there could be better choices than the two suggested ones."
   */

  HYPRE_Real  strong_th = 0.5;  /* 2d=>0.25 (default) 3d=>0.5 */

  HYPRE_BoomerAMGSetStrongThreshold(amg, strong_th);
  HYPRE_BoomerAMGSetStrongThresholdR(amg, strong_th);

  if (amg_as_precond) {

    /* Maximum size of coarsest grid. The default is 9 */

    HYPRE_BoomerAMGSetMaxCoarseSize(amg, 50);

    /* The default smoother is a good compromise in case of preconditioner
     * Remark: fcf-jacobi or l1scaled-jacobi (or chebyshev) as up/down smoothers
     * could be a good choice
     *
     *  0: Jacobi
     *  6: hybrid sym. Gauss Seidel/Jacobi
     *  7: Jacobi
     *  8: l1 sym. Gauss Seidel
     * 13: l1-Gauss Seidel (forward) --> default on the down cycle
     * 14: l1 Gauss Seidel (backward) --> default on the up cycle
     * 15: CG
     * 16: Chebyshev
     * 17: FCF-Jacobi
     * 18: l1 Jacobi
     */

    HYPRE_BoomerAMGSetCycleRelaxType(amg, 13, 1); /* down cycle */
    HYPRE_BoomerAMGSetCycleRelaxType(amg, 14, 2); /* up cycle */
    HYPRE_BoomerAMGSetCycleRelaxType(amg,  8, 3); /* coarsest level */

    /* Number of smoothing steps (precond, num, down/up/coarse) */

    HYPRE_BoomerAMGSetCycleNumSweeps(amg, 1, 1); /* down cycle */
    HYPRE_BoomerAMGSetCycleNumSweeps(amg, 1, 2); /* up cycle */
    HYPRE_BoomerAMGSetCycleNumSweeps(amg, 1, 3); /* coarsest level */

    HYPRE_BoomerAMGSetTol(amg, 0.0);
    HYPRE_BoomerAMGSetMaxIter(amg, 1);

    /* Aggressive coarsening on the nl=2 first levels */

    HYPRE_BoomerAMGSetAggNumLevels(amg, 2);

    /* From HYPRE documentation: For levels with aggressive coarsening
     *
     * AggInterType: interpolation used on levels of aggressive coarsening
     *
     *  4: (default) multipass interpolation
     *  5: 2-stage extended interpolation in matrix-matrix form
     *  6: 2-stage extended+i interpolation in matrix-matrix form
     *  7: 2-stage extended+e interpolation in matrix-matrix form
     */

    HYPRE_BoomerAMGSetAggInterpType(amg, 4);

    /* Number of paths for aggressive coarsening (default = 1) */

    HYPRE_BoomerAMGSetNumPaths(amg, 2);

    /* From the HYPRE documentation.
     * "In order to reduce communication, there is a non-Galerkin coarse grid
     * sparsification option. [...] It is common to drop more aggressively on
     * coarser levels. A common choice of drop-tolerances is [0.0, 0.01, 0.05]
     * where the value of 0.0 will skip the non-Galerkin process on the first
     * coarse level (level 1), use a drop-tolerance of 0.01 on the second coarse
     * level (level 2) and then use 0.05 on all subsequent coarse levels."
     */

    HYPRE_Real  nongalerkin_tol[3] = {0., 0.01, 0.05};
    HYPRE_BoomerAMGSetNonGalerkTol(amg, 3, nongalerkin_tol);

  }
  else { /* AMG as solver */

    HYPRE_BoomerAMGSetRelaxType(amg, 6);      /* Sym G.S./Jacobi hybrid */
    HYPRE_BoomerAMGSetMaxIter(amg, slesp->cvg_param.n_max_iter);
    HYPRE_BoomerAMGSetRelaxOrder(amg, 0);

    /* This option is recommended for GPU usage */

    HYPRE_BoomerAMGSetKeepTranspose(amg, 1);  /* keep transpose to avoid
                                                 SpMTV */
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup hook function for a Hypre solver with preconditioner
 *         This function is called at the end of the setup stage for a solver.
 *
 *         Check HYPRE documentation for available options:
 *         https://hypre.readthedocs.io/en/latest/index.html
 *
 * \param[in]      verbosity   verbosity level
 * \param[in, out] context     pointer to optional (untyped) value or structure
 * \param[in, out] solver_p    pointer to a HYPRE solver (cast ont-the-fly)
 */
/*----------------------------------------------------------------------------*/

static void
_hypre_generic_pc_hook(int     verbosity,
                       void   *context,
                       void   *solver_p)
{
  CS_NO_WARN_IF_UNUSED(verbosity);

  cs_param_sles_t  *slesp = (cs_param_sles_t *)context;
  HYPRE_Solver  hs = solver_p;
  HYPRE_Solver  pc = _set_hypre_solver(slesp, hs);

  switch (slesp->precond) {

  case CS_PARAM_PRECOND_NONE:
  case CS_PARAM_PRECOND_DIAG:
    break;

  case CS_PARAM_PRECOND_ILU0:
    HYPRE_ILUSetMaxIter(pc, (HYPRE_Int)1);
    HYPRE_ILUSetTol(pc, (HYPRE_Real)0.0);
    HYPRE_ILUSetType(pc, (HYPRE_Int)0);
    break;

  case CS_PARAM_PRECOND_BJACOB_ILU0:
    HYPRE_EuclidSetLevel(pc, (HYPRE_Int)0);
    HYPRE_EuclidSetBJ(pc, (HYPRE_Int)1);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: System: %s\n Invalid solver/preconditioner with HYPRE.",
              __func__, slesp->name);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set parameters for initializing SLES structures used for the
 *        resolution of the linear system.
 *        Case of Hypre families of solvers.
 *
 *        Check HYPRE documentation for available options:
 *        https://hypre.readthedocs.io/en/latest/index.html
 *
 * \param[in]      use_field_id   if false use system name
 * \param[in, out] slesp          pointer to a \ref cs_param_sles_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_set_hypre_sles(bool                 use_field_id,
                cs_param_sles_t     *slesp)
{
  assert(slesp != NULL);  /* Sanity checks */

  const char  errmsg[] = "Invalid couple (solver,preconditionner) with HYPRE.";
  const  char  *sles_name = use_field_id ? NULL : slesp->name;
  assert(slesp->field_id > -1 || sles_name != NULL);

  /* Set the solver and the associated preconditioner */

  switch(slesp->solver) {

  case CS_PARAM_ITSOL_AMG:
    cs_sles_hypre_define(slesp->field_id,
                         sles_name,
                         CS_SLES_HYPRE_BOOMERAMG,
                         CS_SLES_HYPRE_NONE,
                         _hypre_boomeramg_hook,
                         (void *)slesp);
    break;

  case CS_PARAM_ITSOL_BICG:
  case CS_PARAM_ITSOL_BICGSTAB2:
    switch (slesp->precond) {

    case CS_PARAM_PRECOND_AMG:
      cs_sles_hypre_define(slesp->field_id,
                           sles_name,
                           CS_SLES_HYPRE_BICGSTAB,
                           CS_SLES_HYPRE_BOOMERAMG,
                           _hypre_boomeramg_hook,
                           (void *)slesp);
      break;

    case CS_PARAM_PRECOND_DIAG:
    case CS_PARAM_PRECOND_NONE:
      cs_sles_hypre_define(slesp->field_id,
                           sles_name,
                           CS_SLES_HYPRE_BICGSTAB,
                           CS_SLES_HYPRE_NONE,
                           _hypre_generic_pc_hook,
                           (void *)slesp);
      break;

    case CS_PARAM_PRECOND_BJACOB_ILU0:
      cs_sles_hypre_define(slesp->field_id,
                           sles_name,
                           CS_SLES_HYPRE_BICGSTAB,
                           CS_SLES_HYPRE_EUCLID,
                           _hypre_generic_pc_hook,
                           (void *)slesp);
      break;

    case CS_PARAM_PRECOND_ILU0:
      cs_sles_hypre_define(slesp->field_id,
                           sles_name,
                           CS_SLES_HYPRE_BICGSTAB,
                           CS_SLES_HYPRE_ILU,
                           _hypre_generic_pc_hook,
                           (void *)slesp);
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " %s: System: %s\n %s\n", __func__, slesp->name, errmsg);
    }
    break;

  case CS_PARAM_ITSOL_CG:
  case CS_PARAM_ITSOL_FCG:
    switch (slesp->precond) {

    case CS_PARAM_PRECOND_AMG:
      cs_sles_hypre_define(slesp->field_id,
                           sles_name,
                           CS_SLES_HYPRE_PCG,
                           CS_SLES_HYPRE_BOOMERAMG,
                           _hypre_boomeramg_hook,
                           (void *)slesp);
      break;

    case CS_PARAM_PRECOND_DIAG:
    case CS_PARAM_PRECOND_NONE:
      cs_sles_hypre_define(slesp->field_id,
                           sles_name,
                           CS_SLES_HYPRE_PCG,
                           CS_SLES_HYPRE_NONE,
                           _hypre_generic_pc_hook,
                           (void *)slesp);
      break;

    case CS_PARAM_PRECOND_BJACOB_ILU0:
      cs_sles_hypre_define(slesp->field_id,
                           sles_name,
                           CS_SLES_HYPRE_PCG,
                           CS_SLES_HYPRE_EUCLID,
                           _hypre_generic_pc_hook,
                           (void *)slesp);
      break;

    case CS_PARAM_PRECOND_ILU0:
      cs_sles_hypre_define(slesp->field_id,
                           sles_name,
                           CS_SLES_HYPRE_PCG,
                           CS_SLES_HYPRE_ILU,
                           _hypre_generic_pc_hook,
                           (void *)slesp);
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " %s: System: %s\n %s\n", __func__, slesp->name, errmsg);
    }
    break;

  case CS_PARAM_ITSOL_FGMRES:
  case CS_PARAM_ITSOL_GCR:
    switch (slesp->precond) {

    case CS_PARAM_PRECOND_AMG:
      cs_sles_hypre_define(slesp->field_id,
                           sles_name,
                           CS_SLES_HYPRE_FLEXGMRES,
                           CS_SLES_HYPRE_BOOMERAMG,
                           _hypre_boomeramg_hook,
                           (void *)slesp);
      break;

    case CS_PARAM_PRECOND_NONE:
    case CS_PARAM_PRECOND_DIAG:
      cs_sles_hypre_define(slesp->field_id,
                           sles_name,
                           CS_SLES_HYPRE_FLEXGMRES,
                           CS_SLES_HYPRE_NONE,
                           _hypre_generic_pc_hook,
                           (void *)slesp);
      break;

    case CS_PARAM_PRECOND_BJACOB_ILU0:
      cs_sles_hypre_define(slesp->field_id,
                           sles_name,
                           CS_SLES_HYPRE_FLEXGMRES,
                           CS_SLES_HYPRE_EUCLID,
                           _hypre_generic_pc_hook,
                           (void *)slesp);
      break;

    case CS_PARAM_PRECOND_ILU0:
      cs_sles_hypre_define(slesp->field_id,
                           sles_name,
                           CS_SLES_HYPRE_FLEXGMRES,
                           CS_SLES_HYPRE_ILU,
                           _hypre_generic_pc_hook,
                           (void *)slesp);
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " %s: System: %s\n Invalid solver/preconditioner with HYPRE.",
                __func__, slesp->name);
    }
    break;

  case CS_PARAM_ITSOL_GMRES:
    switch (slesp->precond) {

    case CS_PARAM_PRECOND_AMG:
      cs_sles_hypre_define(slesp->field_id,
                           sles_name,
                           CS_SLES_HYPRE_GMRES,
                           CS_SLES_HYPRE_BOOMERAMG,
                           _hypre_boomeramg_hook,
                           (void *)slesp);
      break;

    case CS_PARAM_PRECOND_DIAG:
    case CS_PARAM_PRECOND_NONE:
      cs_sles_hypre_define(slesp->field_id,
                           sles_name,
                           CS_SLES_HYPRE_GMRES,
                           CS_SLES_HYPRE_NONE,
                           _hypre_generic_pc_hook,
                           (void *)slesp);
      break;

    case CS_PARAM_PRECOND_BJACOB_ILU0:
      cs_sles_hypre_define(slesp->field_id,
                           sles_name,
                           CS_SLES_HYPRE_GMRES,
                           CS_SLES_HYPRE_EUCLID,
                           _hypre_generic_pc_hook,
                           (void *)slesp);
      break;

    case CS_PARAM_PRECOND_ILU0:
      cs_sles_hypre_define(slesp->field_id,
                           sles_name,
                           CS_SLES_HYPRE_GMRES,
                           CS_SLES_HYPRE_ILU,
                           _hypre_generic_pc_hook,
                           (void *)slesp);
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " %s: System: %s\n Invalid solver/preconditioner with HYPRE.",
                __func__, slesp->name);
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: System: %s\n Incompatible solver with HYPRE.",
              __func__, slesp->name);
    break;

  } /* Switch on solver */
}
#endif /* HAVE_HYPRE */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a new structure storing a set of parameters used when calling
 *        MUMPS. Set a default value for the advanced parameters.
 *
 * \param[in] is_single   single-precision or double-precision
 * \param[in] facto_type  type of factorization to consider
 */
/*----------------------------------------------------------------------------*/

static cs_param_sles_mumps_t *
_create_mumps_param(bool                            is_single,
                    cs_param_sles_facto_type_t      facto_type)
{
  cs_param_sles_mumps_t  *mumpsp = NULL;

  BFT_MALLOC(mumpsp, 1, cs_param_sles_mumps_t);

  mumpsp->is_single = is_single;
  mumpsp->facto_type = facto_type;

  /* Advanced options (default settings) */

  mumpsp->analysis_algo = CS_PARAM_SLES_ANALYSIS_AUTO;
  mumpsp->mem_usage = CS_PARAM_SLES_MEMORY_AUTO;

  mumpsp->advanced_optim = false;   /* No advanced MPI/OpenMP optimization */
  mumpsp->blr_threshold = 0;        /* No BLR */
  mumpsp->mem_coef = -1;            /* No additional memory range */
  mumpsp->block_analysis = 0;       /* No clustered analysis */
  mumpsp->ir_steps = 0;             /* No iterative refinement */

  return mumpsp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy into a new structure the given set of parameters used when
 *        calling MUMPS
 *
 * \param[in] mumpsp   set of mumps parameters
 */
/*----------------------------------------------------------------------------*/

static cs_param_sles_mumps_t *
_copy_mumps_param(const cs_param_sles_mumps_t   *mumpsp)
{
  cs_param_sles_mumps_t  *cpy = NULL;

  BFT_MALLOC(cpy, 1, cs_param_sles_mumps_t);

  cpy->analysis_algo = mumpsp->analysis_algo;
  cpy->facto_type = mumpsp->facto_type;

  cpy->is_single = mumpsp->is_single;
  cpy->mem_usage = mumpsp->mem_usage;
  cpy->advanced_optim = mumpsp->advanced_optim;
  cpy->blr_threshold = mumpsp->blr_threshold;
  cpy->mem_coef = mumpsp->mem_coef;
  cpy->block_analysis = mumpsp->block_analysis;
  cpy->ir_steps = mumpsp->ir_steps;

  return cpy;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log a set of parameters used when calling MUMPS
 *
 * \param[in] name     name related to the current SLES
 * \param[in] mumpsp   set of mumps parameters
 */
/*----------------------------------------------------------------------------*/

static void
_log_mumps_param(const char                    *name,
                 const cs_param_sles_mumps_t   *mumpsp)
{
  char type = (mumpsp->is_single) ? 's' : 'd';
  char tag[32];

  switch (mumpsp->facto_type) {

  case CS_PARAM_SLES_FACTO_LU:
    sprintf(tag, "%cmumps_lu", type);
    break;
  case CS_PARAM_SLES_FACTO_LDLT_SYM:
    sprintf(tag, "%cmumps_ldlt_sym", type);
    break;
  case CS_PARAM_SLES_FACTO_LDLT_SPD:
    sprintf(tag, "%cmumps_ldlt_spd", type);
    break;

  default:
    sprintf(tag, "undefined");
    break;

  }

  cs_log_printf(CS_LOG_SETUP, "  * %s | MUMPS_type:              %s\n",
                name, tag);

  /* Strategy for the memory usage */

  switch (mumpsp->mem_usage) {
  case CS_PARAM_SLES_MEMORY_CONSTRAINED:
    cs_log_printf(CS_LOG_SETUP, "  * %s | Memory_usage:            %s\n",
                  name, "constrained");
    break;
  case CS_PARAM_SLES_MEMORY_AUTO:
    cs_log_printf(CS_LOG_SETUP, "  * %s | Memory_usage:            %s\n",
                  name, "automatic");
    break;
  case CS_PARAM_SLES_MEMORY_CPU_DRIVEN:
    cs_log_printf(CS_LOG_SETUP, "  * %s | Memory_usage:            %s\n",
                  name, "CPU-driven (efficiency first)");
    break;

  default:
    cs_log_printf(CS_LOG_SETUP, "  * %s | Memory_usage:            %s\n",
                  name, "Undefined");
    break;
  }

  /* Algorithm used for the analysis step */

  switch (mumpsp->analysis_algo) {
  case CS_PARAM_SLES_ANALYSIS_AMD:
    cs_log_printf(CS_LOG_SETUP, "  * %s | Analysis_algo:           %s\n",
                  name, "AMD");
    break;
  case CS_PARAM_SLES_ANALYSIS_QAMD:
    cs_log_printf(CS_LOG_SETUP, "  * %s | Analysis_algo:           %s\n",
                  name, "QAMD");
    break;
  case CS_PARAM_SLES_ANALYSIS_PORD:
    cs_log_printf(CS_LOG_SETUP, "  * %s | Analysis_algo:           %s\n",
                  name, "PORD");
    break;
  case CS_PARAM_SLES_ANALYSIS_SCOTCH:
    cs_log_printf(CS_LOG_SETUP, "  * %s | Analysis_algo:           %s\n",
                  name, "SCOTCH");
    break;
  case CS_PARAM_SLES_ANALYSIS_PTSCOTCH:
    cs_log_printf(CS_LOG_SETUP, "  * %s | Analysis_algo:           %s\n",
                  name, "PT-SCOTCH");
    break;
  case CS_PARAM_SLES_ANALYSIS_METIS:
    cs_log_printf(CS_LOG_SETUP, "  * %s | Analysis_algo:           %s\n",
                  name, "METIS");
    break;
  case CS_PARAM_SLES_ANALYSIS_PARMETIS:
    cs_log_printf(CS_LOG_SETUP, "  * %s | Analysis_algo:           %s\n",
                  name, "PARMETIS");
    break;
  case CS_PARAM_SLES_ANALYSIS_AUTO:
    cs_log_printf(CS_LOG_SETUP, "  * %s | Analysis_algo:           %s\n",
                  name, "automatic choice done by MUMPS");
    break;

  default:
    cs_log_printf(CS_LOG_SETUP, "  * %s | Analysis_algo:           %s\n",
                  name, "Undefined");
    break;
  }

  cs_log_printf(CS_LOG_SETUP, "  * %s | Advanced_Optim:          %s\n",
                name, cs_base_strtf(mumpsp->advanced_optim));

  if (mumpsp->block_analysis > 1)
    cs_log_printf(CS_LOG_SETUP, "  * %s | Block_Size in analysis:   %d\n",
                  name, mumpsp->block_analysis);

  if (cs_math_fabs(mumpsp->ir_steps) > 0)
    cs_log_printf(CS_LOG_SETUP, "  * %s | Iterative_Refinement:     %d\n",
                  name, mumpsp->ir_steps);

  if (fabs(mumpsp->blr_threshold) > FLT_MIN)
    cs_log_printf(CS_LOG_SETUP, "  * %s | BLR_threshold:            %e\n",
                  name, mumpsp->blr_threshold);

  if (mumpsp->mem_coef > 0)
    cs_log_printf(CS_LOG_SETUP, "  * %s | Memory pct. increase:    %f\n",
                  name, mumpsp->mem_coef);
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a cs_param_sles_saddle_t structure and assign a minimalist
 *        default settings
 *
 * \return a pointer to the new cs_param_sles_saddle_t structure
 */
/*----------------------------------------------------------------------------*/

cs_param_sles_saddle_t *
cs_param_sles_saddle_create(void)
{
  cs_param_sles_saddle_t  *saddlep = NULL;

  BFT_MALLOC(saddlep, 1, cs_param_sles_saddle_t);

  saddlep->verbosity = 0;
  saddlep->solver = CS_PARAM_SADDLE_SOLVER_NONE;  /* Not used */
  saddlep->precond = CS_PARAM_SADDLE_PRECOND_NONE;
  saddlep->schur_approximation = CS_PARAM_SCHUR_NONE;

  saddlep->cvg_param =  (cs_param_sles_cvg_t) {
    .n_max_iter = 50,
    .atol = 1e-12,       /* absolute tolerance */
    .rtol = 1e-6,        /* relative tolerance */
    .dtol = 1e3 };       /* divergence tolerance */

  saddlep->schur_sles_param = NULL;

  return saddlep;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize a \ref cs_param_sles_t structure for the Schur
 *        approximation nested inside a ref cs_param_sles_saddle_t
 *        structure. By default, this member is not allocated. Do nothing if
 *        the related structure is already allocated.
 *
 * \param[in]      basename   prefix for the naming of the Schur system
 * \param[in, out] saddlep    pointer to the structure to update
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_saddle_init_schur(const char                *basename,
                                cs_param_sles_saddle_t    *saddlep)
{
  if (saddlep == NULL)
    return;

  if (saddlep->schur_sles_param != NULL)
    return; /* Initialization has already been performed */

  char  *schur_name = NULL;
  size_t  len = sizeof(basename) + sizeof("_schur_system");

  BFT_MALLOC(schur_name, len + 1, char);
  sprintf(schur_name, "%s_schur_system", basename);

  cs_param_sles_t  *schur_slesp = cs_param_sles_create(-1, schur_name);

  schur_slesp->precond = CS_PARAM_PRECOND_AMG;   /* preconditioner */
  schur_slesp->solver = CS_PARAM_ITSOL_FCG;      /* iterative solver */
  schur_slesp->amg_type = CS_PARAM_AMG_HOUSE_K;  /* no predefined AMG type */
  schur_slesp->cvg_param.rtol = 1e-4;            /* relative tolerance to stop
                                                    the iterative solver */

  BFT_FREE(schur_name);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy a cs_param_sles_saddle_t structure from ref to dest
 *
 * \param[in]      ref     reference structure to be copied
 * \param[in, out] dest    destination structure
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_saddle_copy(const cs_param_sles_saddle_t   *ref,
                          cs_param_sles_saddle_t         *dest)
{
  if (ref == NULL)
    return;

  dest->solver = ref->solver;
  dest->precond = ref->precond;
  dest->schur_approximation = ref->schur_approximation;

  dest->cvg_param.rtol = ref->cvg_param.rtol;
  dest->cvg_param.atol = ref->cvg_param.atol;
  dest->cvg_param.dtol = ref->cvg_param.dtol;
  dest->cvg_param.n_max_iter = ref->cvg_param.n_max_iter;

  if (ref->schur_sles_param != NULL) {

    if (dest->schur_sles_param == NULL)
      cs_param_sles_saddle_init_schur("automatic", dest);

    cs_param_sles_copy_from(ref->schur_sles_param, dest->schur_sles_param);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free the structure storing the parameter settings for a saddle-point
 *        system
 *
 * \param[in, out] p_saddlep    double pointer to the structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_saddle_free(cs_param_sles_saddle_t    **p_saddlep)
{
  if (p_saddlep == NULL)
    return;

  cs_param_sles_saddle_t  *saddlep = *p_saddlep;

  cs_param_sles_free(&(saddlep->schur_sles_param));

  BFT_FREE(saddlep);
  *p_saddlep = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a \ref cs_param_sles_t structure and assign a default
 *         settings
 *
 * \param[in] field_id      id related to to the variable field or -1
 * \param[in] system_name   name of the system to solve or NULL
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

  slesp->name = NULL;
  if (system_name != NULL) {
    size_t  len = strlen(system_name);
    BFT_MALLOC(slesp->name, len + 1, char);
    strncpy(slesp->name, system_name, len + 1);
  }

  slesp->field_id = field_id;                   /* associated field id */
  slesp->verbosity = 0;                         /* SLES verbosity */
  slesp->setup_done = false;

  slesp->solver_class = CS_PARAM_SLES_CLASS_CS; /* solver family */
  slesp->precond = CS_PARAM_PRECOND_DIAG;       /* preconditioner */
  slesp->solver = CS_PARAM_ITSOL_GCR;           /* iterative solver */
  slesp->flexible = false;                      /* not the flexible variant */
  slesp->restart = 15;                          /* restart after ? iterations */
  slesp->amg_type = CS_PARAM_AMG_NONE;          /* no predefined AMG type */

  slesp->pcd_block_type = CS_PARAM_PRECOND_BLOCK_NONE; /* no block by default */
  slesp->resnorm_type = CS_PARAM_RESNORM_FILTERED_RHS;

  slesp->cvg_param =  (cs_param_sles_cvg_t) {
    .n_max_iter = 10000, /* max. number of iterations */
    .atol = 1e-15,       /* absolute tolerance */
    .rtol = 1e-6,        /* relative tolerance */
    .dtol = 1e3 };       /* divergence tolerance */

  slesp->context_param = NULL;

  return slesp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free a \ref cs_param_sles_t structure
 *
 * \param[in, out] slesp    pointer to a \cs_param_sles_t structure to free
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

  /* One asumes that this context has no pointer to free. This is the case up
     to now, since this process is totally managed by the code. */

  BFT_FREE(slesp->context_param);

  BFT_FREE(slesp);
  slesp = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log information related to the linear settings stored in the
 *        structure
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
    cs_log_printf(CS_LOG_SETUP, "             code_saturne\n");
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

  cs_log_printf(CS_LOG_SETUP, "  * %s | SLES Solver.Name:        %s\n",
                slesp->name, cs_param_get_solver_name(slesp->solver));

  if (slesp->solver == CS_PARAM_ITSOL_MUMPS)
    _log_mumps_param(slesp->name, slesp->context_param);

  else { /* Iterative solvers */

    if (slesp->solver == CS_PARAM_ITSOL_AMG)
      cs_log_printf(CS_LOG_SETUP, "  * %s | SLES AMG.Type:           %s\n",
                    slesp->name, cs_param_get_amg_type_name(slesp->amg_type));

    cs_log_printf(CS_LOG_SETUP, "  * %s | SLES Solver.Precond:     %s\n",
                  slesp->name, cs_param_get_precond_name(slesp->precond));

    if (slesp->precond == CS_PARAM_PRECOND_AMG)
      cs_log_printf(CS_LOG_SETUP, "  * %s | SLES AMG.Type:           %s\n",
                    slesp->name, cs_param_get_amg_type_name(slesp->amg_type));
    else if (slesp->precond == CS_PARAM_PRECOND_MUMPS)
      _log_mumps_param(slesp->name, slesp->context_param);

    cs_log_printf(CS_LOG_SETUP, "  * %s | SLES Block.Precond:      %s\n",
                  slesp->name,
                  cs_param_get_precond_block_name(slesp->pcd_block_type));

    cs_log_printf(CS_LOG_SETUP, "  * %s | SLES Solver.max_iter:    %d\n",
                  slesp->name, slesp->cvg_param.n_max_iter);
    cs_log_printf(CS_LOG_SETUP, "  * %s | SLES Solver.rtol:       % -10.6e\n",
                  slesp->name, slesp->cvg_param.rtol);
    cs_log_printf(CS_LOG_SETUP, "  * %s | SLES Solver.atol:       % -10.6e\n",
                  slesp->name, slesp->cvg_param.atol);

    if (slesp->solver == CS_PARAM_ITSOL_GMRES ||
        slesp->solver == CS_PARAM_ITSOL_FGMRES ||
        slesp->solver == CS_PARAM_ITSOL_GCR)
      cs_log_printf(CS_LOG_SETUP, "  * %s | SLES Solver.Restart:     %d\n",
                    slesp->name, slesp->restart);

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

  } /* Iterative solver */

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
cs_param_sles_copy_from(const cs_param_sles_t   *src,
                        cs_param_sles_t         *dst)
{
  if (src == NULL || dst == NULL)
    return;

  /* Remark: name is managed at the creation of the structure */

  dst->setup_done = src->setup_done;
  dst->verbosity = src->verbosity;
  dst->field_id = src->field_id;

  dst->solver_class = src->solver_class;
  dst->precond = src->precond;
  dst->solver = src->solver;
  dst->amg_type = src->amg_type;
  dst->pcd_block_type = src->pcd_block_type;
  dst->resnorm_type = src->resnorm_type;

  dst->cvg_param.rtol = src->cvg_param.rtol;
  dst->cvg_param.atol = src->cvg_param.atol;
  dst->cvg_param.dtol = src->cvg_param.dtol;
  dst->cvg_param.n_max_iter = src->cvg_param.n_max_iter;

  if (dst->precond == CS_PARAM_PRECOND_MUMPS ||
      dst->solver == CS_PARAM_ITSOL_MUMPS) {

    if (dst->context_param != NULL)
      BFT_FREE(dst->context_param);

    dst->context_param = _copy_mumps_param(src->context_param);

  }
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

  _check_settings(slesp);

  /* use_field_id set to true means that slesp->name is used instead of
     field_id to retrieve the associated SLES structure */

  switch (slesp->solver_class) {

  case CS_PARAM_SLES_CLASS_CS: /* code_saturne's solvers */
    _set_saturne_sles(use_field_id, slesp);
    break;

  case CS_PARAM_SLES_CLASS_MUMPS: /* MUMPS sparse direct solvers */
    _set_mumps_sles(use_field_id, slesp);
    break;

#if defined(HAVE_HYPRE)
  case CS_PARAM_SLES_CLASS_HYPRE: /* HYPRE solvers through PETSc or not */
    _set_hypre_sles(use_field_id, slesp);
    break;

  case CS_PARAM_SLES_CLASS_PETSC: /* PETSc solvers */
    _set_petsc_hypre_sles(use_field_id, slesp);
    break;
#else
  case CS_PARAM_SLES_CLASS_HYPRE: /* HYPRE solvers through PETSc */
  case CS_PARAM_SLES_CLASS_PETSC: /* PETSc solvers */
    _set_petsc_hypre_sles(use_field_id, slesp);
    break;
#endif

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
/*!
 * \brief Set the main memebers of a cs_param_sles_mumps_t structure. This
 *        structure is allocated if needed. Other members are kept to their
 *        values.
 *
 * \param[in, out] slesp         pointer to a cs_param_sles_t structure
 * \param[in]      is_single     single-precision or double-precision
 * \param[in]      facto_type    type of factorization to consider
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_mumps(cs_param_sles_t              *slesp,
                    bool                          is_single,
                    cs_param_sles_facto_type_t    facto_type)
{
  if (slesp == NULL)
    return;

  if (slesp->context_param == NULL)
    slesp->context_param = _create_mumps_param(is_single, facto_type);

  else {

    /* One assumes that the existing context structure is related to MUMPS */

    cs_param_sles_mumps_t  *mumpsp = slesp->context_param;

    mumpsp->is_single = is_single;
    mumpsp->facto_type = facto_type;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the main memebers of a cs_param_sles_mumps_t structure. This
 *        structure is allocated if needed. Other members are kept to their
 *        values. Please refer to the MUMPS user guide for more details about
 *        the following advanced options.
 *
 * \param[in, out] slesp            pointer to a cs_param_sles_t structure
 * \param[in]      analysis_algo    algorithm used for the analysis step
 * \param[in]      block_analysis   > 0: fixed block size; 0: nothing
 * \param[in]      mem_coef         percentage increase in the memory workspace
 * \param[in]      blr_threshold    Accuracy in BLR compression (0: not used)
 * \param[in]      ir_steps         0: No, otherwise number of iterations
 * \param[in]      mem_usage        strategy to adopt for the memory usage
 * \param[in]      advanced_optim   activate advanced optimization (MPI/openMP)
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_mumps_advanced(cs_param_sles_t                *slesp,
                             cs_param_sles_analysis_algo_t   analysis_algo,
                             int                             block_analysis,
                             double                          mem_coef,
                             double                          blr_threshold,
                             int                             ir_steps,
                             cs_param_sles_memory_usage_t    mem_usage,
                             bool                            advanced_optim)
{
  if (slesp == NULL)
    return;

  if (slesp->context_param == NULL)
    slesp->context_param = _create_mumps_param(false, CS_PARAM_SLES_FACTO_LU);

  /* One assumes that the existing context structure is related to MUMPS */

  cs_param_sles_mumps_t  *mumpsp = slesp->context_param;

  mumpsp->analysis_algo = analysis_algo;
  mumpsp->block_analysis = block_analysis;
  mumpsp->mem_coef = mem_coef;
  mumpsp->blr_threshold = blr_threshold;
  mumpsp->ir_steps = CS_MAX(ir_steps, -ir_steps);
  mumpsp->mem_usage = mem_usage;
  mumpsp->advanced_optim = advanced_optim;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the settings associated to a cs_sles_t structure and apply
 *        those defined in the given cs_param_sles_t structure.
 *        This function is used only when a first setup has been performed.
 *
 *        One modifies only some specific options like the max. number of
 *        iterations or the relative tolerance
 *
 * \param[in] use_field_id  if false use a name to retrieve the cs_sles_t struc.
 * \param[in] slesp         pointer to a cs_param_sles_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_update_cvg_settings(bool                     use_field_id,
                                  const cs_param_sles_t   *slesp)
{
  if (slesp == NULL)
    return;

  /* Retrieve the sles structure associated to this sles parameters */

  const  char  *sles_name = use_field_id ? NULL : slesp->name;
  assert(slesp->field_id > -1 || sles_name != NULL);

  cs_sles_t  *sles = cs_sles_find(slesp->field_id, sles_name);

  if (sles == NULL)
    return;

  cs_param_sles_cvg_t  cvgp = slesp->cvg_param;

  switch (slesp->solver_class) {

  case CS_PARAM_SLES_CLASS_CS: /* code_saturne's own solvers */
    {
      switch (slesp->solver) {

      case CS_PARAM_ITSOL_AMG:
        {
          cs_multigrid_t  *mg = cs_sles_get_context(sles);
          assert(mg != NULL);

          cs_multigrid_set_max_cycles(mg, cvgp.n_max_iter);
        }
        break;

      case CS_PARAM_ITSOL_GCR:
      case CS_PARAM_ITSOL_GMRES:
        {
          cs_sles_it_t  *itsol = cs_sles_get_context(sles);
          assert(itsol != NULL);

          cs_sles_it_set_n_max_iter(itsol, cvgp.n_max_iter);
          cs_sles_it_set_restart_interval(itsol, slesp->restart);
        }
        break;

      default:
        {
          cs_sles_it_t  *itsol = cs_sles_get_context(sles);
          assert(itsol != NULL);

          cs_sles_it_set_n_max_iter(itsol, cvgp.n_max_iter);
        }
        break;

      } /* which solver */

    }   /* code_saturne class */
    break;

#if defined(HAVE_PETSC)
  case CS_PARAM_SLES_CLASS_PETSC:
    {
      cs_sles_petsc_t  *petsc_ctx = cs_sles_get_context(sles);
      assert(petsc_ctx);

      cs_sles_petsc_set_cvg_criteria(petsc_ctx,
                                     cvgp.rtol, cvgp.atol, cvgp.dtol,
                                     cvgp.n_max_iter);
    }
    break;
#endif

#if defined(HAVE_HYPRE)
  case CS_PARAM_SLES_CLASS_HYPRE:
    {
      cs_sles_hypre_t  *hypre_ctx = cs_sles_get_context(sles);
      assert(hypre_ctx);

      cs_sles_hypre_set_n_max_iter(hypre_ctx, cvgp.n_max_iter);
    }
    break;
#endif

  default:
    /* CS_PARAM_SLES_CLASS_MUMPS => Nothing to do */
    break;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the related solver class from the amg type
 *
 * \param[in]  amg_type    type of AMG to consider
 *
 * \return the related solver class or CS_PARAM_SLES_CLASS_CS
 */
/*----------------------------------------------------------------------------*/

cs_param_sles_class_t
cs_param_sles_get_class_from_amg(cs_param_amg_type_t   amg_type)
{
  switch (amg_type) {

  case CS_PARAM_AMG_HYPRE_BOOMER_V:
  case CS_PARAM_AMG_HYPRE_BOOMER_W:
    return CS_PARAM_SLES_CLASS_HYPRE;

  case CS_PARAM_AMG_PETSC_GAMG_V:
  case CS_PARAM_AMG_PETSC_GAMG_W:
  case CS_PARAM_AMG_PETSC_PCMG:
    return CS_PARAM_SLES_CLASS_PETSC;

  default:
    return CS_PARAM_SLES_CLASS_CS;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check the availability of a solver library and return the requested
 *        one if this is possible or an alternative or CS_PARAM_SLES_N_CLASSES
 *        if no alternative is available.
 *
 * \param[in]       wanted_class  requested class of solvers
 *
 * \return the available solver class related to the requested class
 */
/*----------------------------------------------------------------------------*/

cs_param_sles_class_t
cs_param_sles_check_class(cs_param_sles_class_t   wanted_class)
{
  switch (wanted_class) {

  case CS_PARAM_SLES_CLASS_CS:  /* No issue */
    return CS_PARAM_SLES_CLASS_CS;

  case CS_PARAM_SLES_CLASS_HYPRE:
    /* ------------------------- */
#if defined(HAVE_HYPRE)
    return CS_PARAM_SLES_CLASS_HYPRE;
#else
#if defined(HAVE_PETSC)
#if defined(PETSC_HAVE_HYPRE)
    return CS_PARAM_SLES_CLASS_HYPRE;
#else
    cs_base_warn(__FILE__, __LINE__);
    bft_printf(" Switch to PETSc library since Hypre is not available");
    return CS_PARAM_SLES_CLASS_PETSC; /* Switch to petsc */
#endif
#else
    return CS_PARAM_SLES_N_CLASSES; /* Neither HYPRE nor PETSc */
#endif
#endif

  case CS_PARAM_SLES_CLASS_PETSC:
    /* ------------------------- */
#if defined(HAVE_PETSC)
    return CS_PARAM_SLES_CLASS_PETSC;
#else
    return CS_PARAM_SLES_N_CLASSES;
#endif

  case CS_PARAM_SLES_CLASS_MUMPS:
    /* ------------------------- */
#if defined(HAVE_MUMPS)
    return CS_PARAM_SLES_CLASS_MUMPS;
#else
#if defined(HAVE_PETSC)
#if defined(PETSC_HAVE_MUMPS)
    cs_base_warn(__FILE__, __LINE__);
    bft_printf(" Switch to PETSc library since MUMPS is not available as"
               " a stand-alone library\n");
    return CS_PARAM_SLES_CLASS_PETSC;
#else
    return CS_PARAM_SLES_N_CLASSES;
#endif  /* PETSC_HAVE_MUMPS */
#else
    return CS_PARAM_SLES_N_CLASSES; /* PETSc without MUMPS  */
#endif  /* HAVE_PETSC */
    return CS_PARAM_SLES_N_CLASSES; /* Neither MUMPS nor PETSc */
#endif

  default:
    return CS_PARAM_SLES_N_CLASSES;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if the setting related to the AMG is consistent with the
 *         solver class.
 *
 * \param[in, out] slesp    pointer to a cs_pparam_sles_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_check_amg(cs_param_sles_t   *slesp)
{
  if (slesp == NULL)
    return;
  if (slesp->precond != CS_PARAM_PRECOND_AMG)
    return;

  switch (slesp->solver_class) {

  case CS_PARAM_SLES_CLASS_PETSC:
#if defined(HAVE_PETSC)
#if defined(PETSC_HAVE_HYPRE) /* PETSC with HYPRE */
    if (slesp->amg_type == CS_PARAM_AMG_HOUSE_V ||
        slesp->amg_type == CS_PARAM_AMG_HOUSE_K)
      slesp->amg_type = CS_PARAM_AMG_PETSC_GAMG_V;
#else
    if (slesp->amg_type == CS_PARAM_AMG_HOUSE_V ||
        slesp->amg_type == CS_PARAM_AMG_HOUSE_K ||
        slesp->amg_type == CS_PARAM_AMG_HYPRE_BOOMER_V)
      slesp->amg_type = CS_PARAM_AMG_PETSC_GAMG_V;
    else if (slesp->amg_type == CS_PARAM_AMG_HYPRE_BOOMER_W)
      slesp->amg_type = CS_PARAM_AMG_PETSC_GAMG_W;
#endif
#else  /* PETSC is not available */
    bft_error(__FILE__, __LINE__, 0,
              " %s(): System \"%s\" PETSc is not available.\n"
              " Please check your installation settings.\n",
              __func__, slesp->name);
#endif
    break;

  case CS_PARAM_SLES_CLASS_HYPRE:
#if defined(HAVE_HYPRE)
    if (slesp->amg_type == CS_PARAM_AMG_HOUSE_V ||
        slesp->amg_type == CS_PARAM_AMG_HOUSE_K ||
        slesp->amg_type == CS_PARAM_AMG_PETSC_PCMG ||
        slesp->amg_type == CS_PARAM_AMG_PETSC_GAMG_V)
      slesp->amg_type = CS_PARAM_AMG_HYPRE_BOOMER_V;
    else if (slesp->amg_type == CS_PARAM_AMG_PETSC_GAMG_W)
      slesp->amg_type = CS_PARAM_AMG_HYPRE_BOOMER_W;
#else
 #if defined(HAVE_PETSC)
  #if defined(PETSC_HAVE_HYPRE) /* PETSC with HYPRE */
    if (slesp->amg_type == CS_PARAM_AMG_HOUSE_V ||
        slesp->amg_type == CS_PARAM_AMG_HOUSE_K ||
        slesp->amg_type == CS_PARAM_AMG_PETSC_PCMG ||
        slesp->amg_type == CS_PARAM_AMG_PETSC_GAMG_V)
      slesp->amg_type = CS_PARAM_AMG_HYPRE_BOOMER_V;
    else if (slesp->amg_type == CS_PARAM_AMG_PETSC_GAMG_W)
      slesp->amg_type = CS_PARAM_AMG_HYPRE_BOOMER_W;
  #else  /* PETSc without HYPRE */
    bft_error(__FILE__, __LINE__, 0,
              " %s(): System \"%s\" HYPRE is not available.\n"
              " Please check your installation settings.\n",
              __func__, slesp->name);

  #endif
 #else  /* Neither HYPRE nor PETSC is available */
    bft_error(__FILE__, __LINE__, 0,
              " %s(): System \"%s\" HYPRE/PETSc is not available.\n"
              " Please check your installation settings.\n",
              __func__, slesp->name);
 #endif
#endif
    break;

  case CS_PARAM_SLES_CLASS_CS:
    if (slesp->amg_type == CS_PARAM_AMG_PETSC_PCMG ||
        slesp->amg_type == CS_PARAM_AMG_PETSC_GAMG_V ||
        slesp->amg_type == CS_PARAM_AMG_PETSC_GAMG_W ||
        slesp->amg_type == CS_PARAM_AMG_HYPRE_BOOMER_V ||
        slesp->amg_type == CS_PARAM_AMG_HYPRE_BOOMER_W)
      slesp->amg_type = CS_PARAM_AMG_HOUSE_K;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s(): System \"%s\" Incompatible setting detected.\n"
              " Please check your installation settings.\n",
              __func__, slesp->name);
    break; /* Nothing to do */
  }
}

#if defined(HAVE_PETSC)
/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the command line option for PETSc
 *
 * \param[in]      use_prefix    need a prefix
 * \param[in]      prefix        optional prefix
 * \param[in]      keyword       command keyword
 * \param[in]      keyval        command value
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_petsc_cmd(bool          use_prefix,
                        const char   *prefix,
                        const char   *keyword,
                        const char   *keyval)
{
  _petsc_cmd(use_prefix, prefix, keyword, keyval);
}
#endif

/*----------------------------------------------------------------------------*/

END_C_DECLS
