/*============================================================================
 * Routines to handle the SLES settings
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

#if defined(HAVE_HYPRE)
#include <HYPRE_krylov.h>
#include <HYPRE_parcsr_ls.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_error.h>
#include <bft_mem.h>

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

#include "cs_param_sles_setup.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

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
 *
 * \param[in] solver  solver to test
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

static inline bool
_system_should_be_sym(cs_param_itsol_type_t  solver)
{
  switch (solver) {

  case CS_PARAM_ITSOL_CG:
  case CS_PARAM_ITSOL_FCG:
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
 * \param[in] prefix        prefix name associated to the current SLES
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
 * \param[in] prefix        prefix name associated to the current SLES
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

static void
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

static void
_petsc_pchypre_hook(const char              *prefix,
                    const cs_param_sles_t   *slesp,
                    bool                     is_symm,
                    PC                       pc)
{
#if defined(PETSC_HAVE_HYPRE)

  CS_UNUSED(is_symm);

  assert(prefix != NULL);
  assert(slesp != NULL);
  assert(slesp->precond == CS_PARAM_PRECOND_AMG);

  cs_param_amg_boomer_t  *bamgp = slesp->context_param;
  assert(bamgp != NULL);

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

  char option[32];

  /* From HYPRE documentation: https://hypre.readthedocs.io/en/lastest
   *
   * for three-dimensional diffusion problems, it is recommended to choose a
   * lower complexity coarsening like HMIS or PMIS (coarsening 10 or 8) and
   * combine it with a distance-two interpolation (interpolation 6 or 7), that
   * is also truncated to 4 or 5 elements per row. Additional reduction in
   * complexity and increased scalability can often be achieved using one or
   * two levels of aggressive coarsening.
   */

  /* Usage as preconditioner induces the two following lines */

  _petsc_cmd(true, prefix, "pc_hypre_boomeramg_max_iter", "1");
  _petsc_cmd(true, prefix, "pc_hypre_boomeramg_tol", "0.0");

  /* Coarsen type: Note that the default coarsening is HMIS in HYPRE */

  switch (bamgp->coarsen_algo) {

  case CS_PARAM_AMG_BOOMER_COARSEN_FALGOUT:
    sprintf(option, "%s", "Falgout");
    break;
  case CS_PARAM_AMG_BOOMER_COARSEN_PMIS:
    sprintf(option, "%s", "PMIS");
    break;
  case CS_PARAM_AMG_BOOMER_COARSEN_HMIS:
    sprintf(option, "%s", "HMIS");
    break;

  case CS_PARAM_AMG_BOOMER_COARSEN_CGC:
  case CS_PARAM_AMG_BOOMER_COARSEN_CGC_E:
    bft_error(__FILE__, __LINE__, 0, "%s: Not available from PETSc.", __func__);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Undefined coarsening algo.", __func__);
    break;
  }

  _petsc_cmd(true, prefix, "pc_hypre_boomeramg_coarsen_type", option);

  /* Note that the default interpolation is extended+i interpolation truncated
   * to 4 elements per row. Using 0 means there is no limitation.
   * good choices are: ext+i-cc, ext+i, FF1
   */

  switch (bamgp->interp_algo) {

  case CS_PARAM_AMG_BOOMER_INTERP_HYPERBOLIC:
    bft_error(__FILE__, __LINE__, 0, "%s: Not available from PETSc.", __func__);
    break;
  case CS_PARAM_AMG_BOOMER_INTERP_EXT_PLUS_I_CC:
    sprintf(option, "%s", "ext+i");  /* Strange but the file
                                        cf. src/ksp/pc/impls/hypre/hypre.c */
    break;
  case CS_PARAM_AMG_BOOMER_INTERP_EXT_PLUS_I:
    sprintf(option, "%s", "ext+i-cc"); /* Strange but the file
                                          cf. src/ksp/pc/impls/hypre/hypre.c */
    break;
  case CS_PARAM_AMG_BOOMER_INTERP_FF1:
    sprintf(option, "%s", "FF1");
    break;
  case CS_PARAM_AMG_BOOMER_INTERP_EXTENDED:
    sprintf(option, "%s", "ext");
    break;
  case CS_PARAM_AMG_BOOMER_INTERP_EXT_PLUS_I_MATRIX:
    sprintf(option, "%s", "ext+i-mm");
    break;
  case CS_PARAM_AMG_BOOMER_INTERP_EXT_PLUS_E_MATRIX:
    sprintf(option, "%s", "ext+e-mm");
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Undefined interpol. algo.", __func__);
    break;
  }

  _petsc_cmd(true, prefix, "pc_hypre_boomeramg_interp_type", option);

  sprintf(option, "%d", bamgp->p_max);
  _petsc_cmd(true, prefix, "pc_hypre_boomeramg_P_max", option);

  /* Number of levels (starting from the finest one) on which one applies an
     aggressive coarsening */

  sprintf(option, "%d", bamgp->n_agg_levels);
  _petsc_cmd(true, prefix, "pc_hypre_boomeramg_agg_nl", option);

  /* Number of paths for aggressive coarsening (default = 1) */

  sprintf(option, "%d", bamgp->n_agg_paths);
  _petsc_cmd(true, prefix, "pc_hypre_boomeramg_agg_num_paths", option);

  /* For best performance, it might be necessary to set certain parameters,
   * which will affect both coarsening and interpolation. One important
   * parameter is the strong threshold.  The default value is 0.25, which
   * appears to be a good choice for 2-dimensional problems and the low
   * complexity coarsening algorithms. For 3-dimensional problems a better
   * choice appears to be 0.5, when using the default coarsening
   * algorithm. However, the choice of the strength threshold is problem
   * dependent.
   */

  sprintf(option, "%.3f", bamgp->strong_threshold);
  _petsc_cmd(true, prefix, "pc_hypre_boomeramg_strong_threshold", option);

  _petsc_cmd(true, prefix, "pc_hypre_boomeramg_no_CF","");

  /* _petsc_cmd(true, prefix, "pc_hypre_boomeramg_grid_sweeps_down","2"); */
  /* _petsc_cmd(true, prefix, "pc_hypre_boomeramg_grid_sweeps_up","2"); */
  /* _petsc_cmd(true, prefix, "pc_hypre_boomeramg_smooth_type","Euclid"); */

  /* Remark: FCF-Jacobi or l1scaled-Jacobi (or Chebyshev) as up/down smoothers
   *  can be a good choice
   */

  /* Smoother for the down cycle */

  switch (bamgp->down_smoother) {

  case CS_PARAM_AMG_BOOMER_JACOBI:
    sprintf(option, "%s", "Jacobi");
    break;
  case CS_PARAM_AMG_BOOMER_FORWARD_GS:
    sprintf(option, "%s", "SOR/Jacobi");
    break;
  case CS_PARAM_AMG_BOOMER_BACKWARD_GS:
    sprintf(option, "%s", "backward-SOR/Jacobi");
    break;
  case CS_PARAM_AMG_BOOMER_HYBRID_SSOR:
    sprintf(option, "%s", "symmetric-SOR/Jacobi");
    break;
  case CS_PARAM_AMG_BOOMER_L1_SGS:
    sprintf(option, "%s", "l1scaled-SOR/Jacobi");
    break;
  case CS_PARAM_AMG_BOOMER_GAUSS_ELIM:
    sprintf(option, "%s", "Gaussian-elimination");
    break;
  case CS_PARAM_AMG_BOOMER_BACKWARD_L1_GS:
    sprintf(option, "%s", "l1-Gauss-Seidel");
    break;
  case CS_PARAM_AMG_BOOMER_FORWARD_L1_GS:
    sprintf(option, "%s", "backward-l1-Gauss-Seidel");
    break;
  case CS_PARAM_AMG_BOOMER_CG:
    sprintf(option, "%s", "CG");
    break;
  case CS_PARAM_AMG_BOOMER_CHEBYSHEV:
    sprintf(option, "%s", "Chebyshev");
    break;
  case CS_PARAM_AMG_BOOMER_FCF_JACOBI:
    sprintf(option, "%s", "FCF-Jacobi");
    break;
  case CS_PARAM_AMG_BOOMER_L1_JACOBI:
    sprintf(option, "%s", "l1scaled-Jacobi");
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid down smoother", __func__);
  }

  _petsc_cmd(true, prefix, "pc_hypre_boomeramg_relax_type_down", option);

  sprintf(option, "%d", bamgp->n_down_iter);
  _petsc_cmd(true, prefix, "pc_hypre_boomeramg_grid_sweeps_down", option);

  /* Smoother for the Up cycle */

  switch (bamgp->up_smoother) {

  case CS_PARAM_AMG_BOOMER_JACOBI:
    sprintf(option, "%s", "Jacobi");
    break;
  case CS_PARAM_AMG_BOOMER_FORWARD_GS:
    sprintf(option, "%s", "SOR/Jacobi");
    break;
  case CS_PARAM_AMG_BOOMER_BACKWARD_GS:
    sprintf(option, "%s", "backward-SOR/Jacobi");
    break;
  case CS_PARAM_AMG_BOOMER_HYBRID_SSOR:
    sprintf(option, "%s", "symmetric-SOR/Jacobi");
    break;
  case CS_PARAM_AMG_BOOMER_L1_SGS:
    sprintf(option, "%s", "l1scaled-SOR/Jacobi");
    break;
  case CS_PARAM_AMG_BOOMER_GAUSS_ELIM:
    sprintf(option, "%s", "Gaussian-elimination");
    break;
  case CS_PARAM_AMG_BOOMER_BACKWARD_L1_GS:
    sprintf(option, "%s", "l1-Gauss-Seidel");
    break;
  case CS_PARAM_AMG_BOOMER_FORWARD_L1_GS:
    sprintf(option, "%s", "backward-l1-Gauss-Seidel");
    break;
  case CS_PARAM_AMG_BOOMER_CG:
    sprintf(option, "%s", "CG");
    break;
  case CS_PARAM_AMG_BOOMER_CHEBYSHEV:
    sprintf(option, "%s", "Chebyshev");
    break;
  case CS_PARAM_AMG_BOOMER_FCF_JACOBI:
    sprintf(option, "%s", "FCF-Jacobi");
    break;
  case CS_PARAM_AMG_BOOMER_L1_JACOBI:
    sprintf(option, "%s", "l1scaled-Jacobi");
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid up smoother", __func__);
  }

  _petsc_cmd(true, prefix, "pc_hypre_boomeramg_relax_type_up", option);

  sprintf(option, "%d", bamgp->n_up_iter);
  _petsc_cmd(true, prefix, "pc_hypre_boomeramg_grid_sweeps_up", option);

  /* Solver for the coarsest level */

  switch (bamgp->coarse_solver) {

  case CS_PARAM_AMG_BOOMER_JACOBI:
    sprintf(option, "%s", "Jacobi");
    break;
  case CS_PARAM_AMG_BOOMER_FORWARD_GS:
    sprintf(option, "%s", "SOR/Jacobi");
    break;
  case CS_PARAM_AMG_BOOMER_BACKWARD_GS:
    sprintf(option, "%s", "backward-SOR/Jacobi");
    break;
  case CS_PARAM_AMG_BOOMER_HYBRID_SSOR:
    sprintf(option, "%s", "symmetric-SOR/Jacobi");
    break;
  case CS_PARAM_AMG_BOOMER_L1_SGS:
    sprintf(option, "%s", "l1scaled-SOR/Jacobi");
    break;
  case CS_PARAM_AMG_BOOMER_GAUSS_ELIM:
    sprintf(option, "%s", "Gaussian-elimination");
    break;
  case CS_PARAM_AMG_BOOMER_BACKWARD_L1_GS:
    sprintf(option, "%s", "l1-Gauss-Seidel");
    break;
  case CS_PARAM_AMG_BOOMER_FORWARD_L1_GS:
    sprintf(option, "%s", "backward-l1-Gauss-Seidel");
    break;
  case CS_PARAM_AMG_BOOMER_CG:
    sprintf(option, "%s", "CG");
    break;
  case CS_PARAM_AMG_BOOMER_CHEBYSHEV:
    sprintf(option, "%s", "Chebyshev");
    break;
  case CS_PARAM_AMG_BOOMER_FCF_JACOBI:
    sprintf(option, "%s", "FCF-Jacobi");
    break;
  case CS_PARAM_AMG_BOOMER_L1_JACOBI:
    sprintf(option, "%s", "l1scaled-Jacobi");
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid up smoother", __func__);
  }

#endif /* defined(PETSC_HAVE_HYPRE) */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the command line options for PC according to the kind of
 *        preconditioner to apply
 *
 * \param[in, out] slesp  set of parameters for the linear algebra
 * \param[in, out] ksp    PETSc solver structure
 */
/*----------------------------------------------------------------------------*/

static void
_petsc_set_pc_type(cs_param_sles_t  *slesp,
                   KSP               ksp)
{
  if (slesp->solver == CS_PARAM_ITSOL_MUMPS)
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
    if (slesp->solver_class == CS_PARAM_SOLVER_CLASS_HYPRE) {
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
      cs_log_printf(CS_LOG_WARNINGS,
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
      cs_log_printf(CS_LOG_WARNINGS,
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
    if (slesp->solver_class == CS_PARAM_SOLVER_CLASS_HYPRE) {

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
        cs_log_printf(CS_LOG_WARNINGS,
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
    cs_param_sles_setup_petsc_pc_amg(slesp->name, slesp, pc);
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
 * \param[in]      slesp  pointer to SLES parameters
 * \param[in, out] ksp    pointer to PETSc KSP context
 */
/*----------------------------------------------------------------------------*/

static void
_petsc_set_krylov_solver(cs_param_sles_t  *slesp,
                         KSP               ksp)
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

  case CS_PARAM_ITSOL_FCG:       /* Flexible Conjugate Gradient */
    KSPSetType(ksp, KSPFCG);
    break;

  case CS_PARAM_ITSOL_FGMRES:    /* Preconditioned flexible GMRES */
    KSPSetType(ksp, KSPFGMRES);
    break;

  case CS_PARAM_ITSOL_GCR:       /* Generalized Conjugate Residual */
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

  case CS_PARAM_ITSOL_GMRES:  /* Preconditioned GMRES */
  case CS_PARAM_ITSOL_FGMRES: /* Flexible GMRES */
    KSPGMRESSetRestart(ksp, slesp->restart);
    break;

  case CS_PARAM_ITSOL_GCR: /* Preconditioned GCR */
    KSPGCRSetRestart(ksp, slesp->restart);
    break;

  case CS_PARAM_ITSOL_MUMPS:
#if defined(PETSC_HAVE_MUMPS)
    {
      cs_param_mumps_t  *mumpsp = slesp->context_param;
      assert(mumpsp != NULL);

      PC  pc;
      KSPGetPC(ksp, &pc);

      if (mumpsp->facto_type == CS_PARAM_MUMPS_FACTO_LU) {

        PCSetType(pc, PCLU);
        PCFactorSetMatSolverType(pc, MATSOLVERMUMPS);

      }
      else {

        assert(mumpsp->facto_type == CS_PARAM_MUMPS_FACTO_LDLT_SPD ||
               mumpsp->facto_type == CS_PARAM_MUMPS_FACTO_LDLT_SYM);

        if (mumpsp->facto_type == CS_PARAM_MUMPS_FACTO_LDLT_SPD) {

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
#else
    bft_error(__FILE__, __LINE__, 0,
              "%s: MUMPS inside PETSc is not available.\n",
              " Please check your settings or your installation.", __func__);
#endif
    break;


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
_petsc_setup_hook(void  *context,
                  void  *ksp_struct)
{
  cs_param_sles_t  *slesp = (cs_param_sles_t  *)context;
  KSP  ksp = ksp_struct;

  cs_fp_exception_disable_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */

  cs_param_sles_setup_petsc_ksp(slesp->name, slesp, ksp);

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

  switch (slesp->precond_block_type) {
  case CS_PARAM_PRECOND_BLOCK_UPPER_TRIANGULAR:
  case CS_PARAM_PRECOND_BLOCK_LOWER_TRIANGULAR:
    PCFieldSplitSetType(pc, PC_COMPOSITE_MULTIPLICATIVE);
    break;

  case CS_PARAM_PRECOND_BLOCK_SYM_GAUSS_SEIDEL:
    PCFieldSplitSetType(pc, PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE);
    break;

  case CS_PARAM_PRECOND_BLOCK_DIAG:
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
 * \brief Function pointer: setup hook for setting PETSc solver and
 *        preconditioner.
 *        Case of multiplicative AMG block preconditioner for a CG with GAMG
 *        as AMG type
 *
 * \param[in, out] context    pointer to optional (untyped) value or structure
 * \param[in, out] ksp_struct pointer to PETSc KSP context
 */
/*----------------------------------------------------------------------------*/

static void
_petsc_amg_block_gamg_hook(void  *context,
                           void  *ksp_struct)
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

  cs_fp_exception_restore_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function pointer: setup hook for setting PETSc solver and
 *        preconditioner.
 *        Case of multiplicative AMG block preconditioner for a CG with boomer
 *        as AMG type
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

  PCSetFromOptions(pc);
  KSPSetFromOptions(ksp);

  cs_fp_exception_restore_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */
}

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
      if (slesp->solver_class == CS_PARAM_SOLVER_CLASS_HYPRE) {
        if (cs_param_sles_hypre_from_petsc()) {

          _petsc_cmd(true, prefix, "ksp_type", "preonly");
          _petsc_cmd(true, prefix, "pc_type", "hypre");
          _petsc_cmd(true, prefix, "pc_hypre_type", "euclid");
          _petsc_cmd(true, prefix, "pc_hypre_euclid_level", "0");

        }
        else
          bft_error(__FILE__, __LINE__, 0,
                    " %s: Invalid option: HYPRE is not installed.", __func__);
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
    case CS_PARAM_PRECOND_MUMPS:
      _petsc_cmd(true, prefix, "ksp_type", "preonly");
#if defined(PETSC_HAVE_MUMPS)
      {
        assert(slesp != NULL);
        assert(slesp->context_param != NULL);

        cs_param_mumps_t  *mumpsp = slesp->context_param;

        if (mumpsp->facto_type == CS_PARAM_MUMPS_FACTO_LDLT_SPD)
          _petsc_cmd(true, prefix, "pc_type", "cholesky");
        else
          _petsc_cmd(true, prefix, "pc_type", "lu");

        _petsc_cmd(true, prefix, "pc_factor_mat_solver_type", "mumps");

        switch (mumpsp->analysis_algo) {

        case CS_PARAM_MUMPS_ANALYSIS_AMD:
          _petsc_cmd(true, prefix, "mat_mumps_icntl_28", "1");
          _petsc_cmd(true, prefix, "mat_mumps_icntl_7", "0");
          break;

        case CS_PARAM_MUMPS_ANALYSIS_QAMD:
          _petsc_cmd(true, prefix, "mat_mumps_icntl_28", "1");
          _petsc_cmd(true, prefix, "mat_mumps_icntl_7", "6");
          break;

        case CS_PARAM_MUMPS_ANALYSIS_PORD:
          _petsc_cmd(true, prefix, "mat_mumps_icntl_28", "1");
          _petsc_cmd(true, prefix, "mat_mumps_icntl_7", "4");
          break;

        case CS_PARAM_MUMPS_ANALYSIS_SCOTCH:
          _petsc_cmd(true, prefix, "mat_mumps_icntl_28", "1");
          _petsc_cmd(true, prefix, "mat_mumps_icntl_7", "3");
          _petsc_cmd(true, prefix, "mat_mumps_icntl_58", "2");
          break;

        case CS_PARAM_MUMPS_ANALYSIS_PTSCOTCH:
          _petsc_cmd(true, prefix, "mat_mumps_icntl_28", "2");
          _petsc_cmd(true, prefix, "mat_mumps_icntl_29", "1");
          _petsc_cmd(true, prefix, "mat_mumps_icntl_58", "0");
          break;

        case CS_PARAM_MUMPS_ANALYSIS_METIS:
          _petsc_cmd(true, prefix, "mat_mumps_icntl_28", "1");
          _petsc_cmd(true, prefix, "mat_mumps_icntl_7", "5");
          break;

        case CS_PARAM_MUMPS_ANALYSIS_PARMETIS:
          _petsc_cmd(true, prefix, "mat_mumps_icntl_28", "2");
          _petsc_cmd(true, prefix, "mat_mumps_icntl_29", "2");
          _petsc_cmd(true, prefix, "mat_mumps_icntl_58", "2");
          break;

        default: /* CS_PARAM_MUMPS_ANALYSIS_AUTO: */
          _petsc_cmd(true, prefix, "mat_mumps_icntl_7", "7");
          break;

        } /* Type of algorithm for the analysis step */

      }
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

  PCSetFromOptions(pc);
  KSPSetFromOptions(ksp);

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
  if (slesp->solver == CS_PARAM_ITSOL_MUMPS) { /* Checks related to MUMPS */

    cs_param_solver_class_t  ret_class =
      cs_param_sles_check_class(CS_PARAM_SOLVER_CLASS_MUMPS);
    if (ret_class == CS_PARAM_N_SOLVER_CLASSES)
      bft_error(__FILE__, __LINE__, 0,
                " %s: Error detected while setting the SLES \"%s\"\n"
                " MUMPS is not available with your installation.\n"
                " Please check your installation settings.\n",
                __func__, slesp->name);
    else
      slesp->solver_class = ret_class;

  }
  else {
    if (slesp->solver_class == CS_PARAM_SOLVER_CLASS_MUMPS)
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
 * \param[in]      use_field_id  if false use system name
 * \param[in, out] slesp         pointer to a \ref cs_param_sles_t structure
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

  cs_sles_t  *sles = cs_sles_find_or_add(slesp->field_id, sles_name);
  cs_sles_it_t  *itsol = cs_sles_get_context(sles);

  if (itsol != NULL && slesp->field_id > -1)
    return; /* Solver settings already forced */

  int  poly_degree = _get_poly_degree(slesp);

  /* Retrieve associated context structures */

  cs_sles_pc_t  *pc = cs_sles_it_get_pc(itsol);
  cs_multigrid_t  *mg = NULL;

  /* 1- Define the iterative solver
   *    =========================== */

  if (itsol == NULL) { /* Not already defined. Add a new one. */

    switch (slesp->solver) {

    case CS_PARAM_ITSOL_AMG:
      {
        switch (slesp->amg_type) {

        case CS_PARAM_AMG_INHOUSE_V:
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

        case CS_PARAM_AMG_INHOUSE_K:
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

    case CS_PARAM_ITSOL_GAUSS_SEIDEL:
      itsol = cs_sles_it_define(slesp->field_id, sles_name,
                                CS_SLES_P_GAUSS_SEIDEL,
                                -1, /* Not useful to apply a preconditioner */
                                slesp->cvg_param.n_max_iter);
      break;

    case CS_PARAM_ITSOL_FGMRES:  /* Not available --> close to GCR */
      cs_base_warn(__FILE__, __LINE__);
      cs_log_printf(CS_LOG_WARNINGS,
                    "%s: Switch to the GCR implementation of code_saturne\n",
                    __func__);
      /* No break (wanted behavior) */
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
                " %s: Invalid iterative solver for solving equation \"%s\".\n"
                " Please modify your settings.", __func__, slesp->name);
      break;

    } /* End of switch */

  }
  else {

    /* The itsol structure has already been defined. Check that this is
       consistent with the SLES parameters */

    cs_sles_it_type_t  itsol_type = cs_sles_it_get_type(itsol);
    int  ret_code = -1;
    switch (itsol_type) {

    case CS_SLES_PCG:
      if (slesp->solver != CS_PARAM_ITSOL_CG)
        ret_code = 0;
      break;

    case CS_SLES_FCG:
    case CS_SLES_IPCG:
      if (slesp->solver != CS_PARAM_ITSOL_FCG)
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
      if (slesp->solver != CS_PARAM_ITSOL_GCR)
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
      cs_log_printf(CS_LOG_WARNINGS,
                    "%s: A flexible Krylov method should be used.\n", __func__);
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
      if (slesp->amg_type == CS_PARAM_AMG_INHOUSE_V) {

        pc = cs_multigrid_pc_create(CS_MULTIGRID_V_CYCLE);
        mg = cs_sles_pc_get_context(pc);
        cs_sles_it_transfer_pc(itsol, &pc);

      }
      else if (slesp->amg_type == CS_PARAM_AMG_INHOUSE_K) {

        pc = cs_multigrid_pc_create(CS_MULTIGRID_K_CYCLE);
        mg = cs_sles_pc_get_context(pc);
        cs_sles_it_transfer_pc(itsol, &pc);

        /* Change the default settings when used as preconditioner */

        cs_multigrid_set_solver_options
          (mg,
           CS_SLES_P_SYM_GAUSS_SEIDEL, /* descent smoother */
           CS_SLES_P_SYM_GAUSS_SEIDEL, /* ascent smoother */
           CS_SLES_PCG,                /* coarse solver */
           1,                          /* n_max_cycles */
           1,                          /* n_max_iter_descent, */
           4,                          /* n_max_iter_ascent */
           500,                        /* n_max_iter_coarse */
           0,                          /* poly_degree_descent */
           0,                          /* poly_degree_ascent */
           1,                          /* poly_degree_coarse */
           -1.0,                       /* precision_mult_descent */
           -1.0,                       /* precision_mult_ascent */
           1.0);                       /* precision_mult_coarse */

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
                            CS_PARAM_MUMPS_FACTO_LU);

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

  /* In case of high verbosity, additional outputs are generated */

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
    cs_param_sles_mumps_reset(slesp);

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

  if (slesp->precond_block_type != CS_PARAM_PRECOND_BLOCK_NONE) {

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

        if (cs_param_sles_hypre_from_petsc())
          cs_sles_petsc_define(slesp->field_id,
                               sles_name,
                               MATMPIAIJ,
                               _petsc_amg_block_boomer_hook,
                               (void *)slesp);

        else {

          cs_base_warn(__FILE__, __LINE__);
          cs_log_printf(CS_LOG_WARNINGS,
                        " %s: System: %s.\n"
                        " Boomer is not available. Switch to GAMG solver.",
                        __func__, slesp->name);

          cs_sles_petsc_define(slesp->field_id,
                               sles_name,
                               MATMPIAIJ,
                               _petsc_amg_block_gamg_hook,
                               (void *)slesp);
        }

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

#if defined(PETSC_HAVE_HYPRE)
    if (slesp->precond == CS_PARAM_PRECOND_AMG) {
      if (slesp->amg_type == CS_PARAM_AMG_HYPRE_BOOMER_V ||
          slesp->amg_type == CS_PARAM_AMG_HYPRE_BOOMER_W) {

        cs_param_amg_boomer_t  *bamgp = slesp->context_param;

        if (bamgp == NULL) /* Define a default set of parameters */
          cs_param_sles_boomeramg_reset(slesp);

      }
    }
#endif

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
  cs_param_amg_boomer_t  *bamgp = slesp->context_param;
  HYPRE_Solver  hs = solver_p;
  HYPRE_Solver  amg = _set_hypre_solver(slesp, hs);
  bool  amg_as_precond = (slesp->solver == CS_PARAM_ITSOL_AMG) ? false : true;

  assert(bamgp != NULL);

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

  HYPRE_BoomerAMGSetCoarsenType(amg, bamgp->coarsen_algo);

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

  HYPRE_BoomerAMGSetInterpType(amg, bamgp->interp_algo);
  HYPRE_BoomerAMGSetPMaxElmts(amg, bamgp->p_max);

  /* From the HYPRE documentation.
   * "One important parameter is the strong threshold, which can be set using
   * the function HYPRE_BoomerAMGSetStrongThreshold. The default value is 0.25,
   * which appears to be a good choice for 2- dimensional problems and the low
   * complexity coarsening algorithms. For 3-dimensional problems a better
   * choice appears to be 0.5, when using the default coarsening algorithm.
   * However, the choice of the strength threshold is problem dependent and
   * therefore there could be better choices than the two suggested ones."
   */

  HYPRE_Real  strong_th = bamgp->strong_threshold;  /* 2d=>0.25 (default)
                                                       3d=>0.5 */

  HYPRE_BoomerAMGSetStrongThreshold(amg, strong_th);
  HYPRE_BoomerAMGSetStrongThresholdR(amg, strong_th);

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

  /* Down cycle: smoother and number of sweeps */

  HYPRE_BoomerAMGSetCycleRelaxType(amg, bamgp->down_smoother, 1);
  HYPRE_BoomerAMGSetCycleNumSweeps(amg, bamgp->n_down_iter, 1);

  /* Up cycle: smoother and number of sweeps */

  HYPRE_BoomerAMGSetCycleRelaxType(amg, bamgp->up_smoother, 2);
  HYPRE_BoomerAMGSetCycleNumSweeps(amg, bamgp->n_up_iter, 2);

  /* Coarsest level */

  HYPRE_BoomerAMGSetCycleRelaxType(amg, bamgp->coarse_solver, 3);

  /* Aggressive coarsening on the nl=2 first levels */

  HYPRE_BoomerAMGSetAggNumLevels(amg, bamgp->n_agg_levels);

  /* Number of paths for aggressive coarsening (default = 1) */

  HYPRE_BoomerAMGSetNumPaths(amg, bamgp->n_agg_paths);

  if (amg_as_precond) {

    /* Maximum size of coarsest grid. The default is 9 */

    HYPRE_BoomerAMGSetMaxCoarseSize(amg, 50);

    HYPRE_BoomerAMGSetTol(amg, 0.0);
    HYPRE_BoomerAMGSetMaxIter(amg, 1);

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

    HYPRE_BoomerAMGSetMaxIter(amg, slesp->cvg_param.n_max_iter);

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
    {
      cs_param_amg_boomer_t  *bamgp = slesp->context_param;

      if (bamgp == NULL) /* Define a default set of parameters */
        cs_param_sles_boomeramg_reset(slesp);

      cs_sles_hypre_define(slesp->field_id,
                           sles_name,
                           CS_SLES_HYPRE_BOOMERAMG,
                           CS_SLES_HYPRE_NONE,
                           _hypre_boomeramg_hook,
                           (void *)slesp);
    }
    break;

  case CS_PARAM_ITSOL_BICG:
  case CS_PARAM_ITSOL_BICGSTAB2:
    switch (slesp->precond) {

    case CS_PARAM_PRECOND_AMG:
      {
        cs_param_amg_boomer_t  *bamgp = slesp->context_param;

        if (bamgp == NULL) /* Define a default set of parameters */
          cs_param_sles_boomeramg_reset(slesp);

        cs_sles_hypre_define(slesp->field_id,
                             sles_name,
                             CS_SLES_HYPRE_BICGSTAB,
                             CS_SLES_HYPRE_BOOMERAMG,
                             _hypre_boomeramg_hook,
                             (void *)slesp);
      }
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
      {
        cs_param_amg_boomer_t  *bamgp = slesp->context_param;

        if (bamgp == NULL) /* Define a default set of parameters */
          cs_param_sles_boomeramg_reset(slesp);

        cs_sles_hypre_define(slesp->field_id,
                             sles_name,
                             CS_SLES_HYPRE_FLEXGMRES,
                             CS_SLES_HYPRE_BOOMERAMG,
                             _hypre_boomeramg_hook,
                             (void *)slesp);
      }
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
      {
        cs_param_amg_boomer_t  *bamgp = slesp->context_param;

        if (bamgp == NULL) /* Define a default set of parameters */
          cs_param_sles_boomeramg_reset(slesp);

        cs_sles_hypre_define(slesp->field_id,
                             sles_name,
                             CS_SLES_HYPRE_GMRES,
                             CS_SLES_HYPRE_BOOMERAMG,
                             _hypre_boomeramg_hook,
                             (void *)slesp);
      }
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

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a cs_sles_t structure in accordance with the settings in a
 *        cs_param_sles_t structure (SLES = Sparse Linear Equation Solver)
 *
 * \param[in]      use_field_id  if false use system name to define a SLES
 * \param[in, out] slesp         pointer to a cs_param_sles_t structure
 *
 * \return an error code (-1 if a problem is encountered, 0 otherwise)
 */
/*----------------------------------------------------------------------------*/

int
cs_param_sles_setup(bool              use_field_id,
                    cs_param_sles_t  *slesp)
{
  if (slesp == NULL)
    return 0;

  _check_settings(slesp);

  /* use_field_id set to true means that slesp->name is used instead of
     field_id to retrieve the associated SLES structure */

  switch (slesp->solver_class) {

  case CS_PARAM_SOLVER_CLASS_CS: /* code_saturne's solvers */
    _set_saturne_sles(use_field_id, slesp);
    break;

  case CS_PARAM_SOLVER_CLASS_MUMPS: /* MUMPS sparse direct solvers */
    _set_mumps_sles(use_field_id, slesp);
    break;

#if defined(HAVE_HYPRE)
  case CS_PARAM_SOLVER_CLASS_HYPRE: /* HYPRE solvers through PETSc or not */
    _set_hypre_sles(use_field_id, slesp);
    break;

  case CS_PARAM_SOLVER_CLASS_PETSC: /* PETSc solvers */
    _set_petsc_hypre_sles(use_field_id, slesp);
    break;
#else
  case CS_PARAM_SOLVER_CLASS_HYPRE: /* HYPRE solvers through PETSc */
  case CS_PARAM_SOLVER_CLASS_PETSC: /* PETSc solvers */
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
cs_param_sles_setup_cvg_param(bool                    use_field_id,
                              const cs_param_sles_t  *slesp)
{
  if (slesp == NULL)
    return;

  /* Retrieve the sles structure associated to this sles parameters */

  const  char  *sles_name = use_field_id ? NULL : slesp->name;
  assert(slesp->field_id > -1 || sles_name != NULL);

  cs_sles_t  *sles = cs_sles_find(slesp->field_id, sles_name);

  if (sles == NULL)
    return;

  cs_param_convergence_t  cvgp = slesp->cvg_param;

  switch (slesp->solver_class) {

  case CS_PARAM_SOLVER_CLASS_CS: /* code_saturne's own solvers */
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

    } /* code_saturne class */
    break;

#if defined(HAVE_PETSC)
  case CS_PARAM_SOLVER_CLASS_PETSC:
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
  case CS_PARAM_SOLVER_CLASS_HYPRE:
    {
      cs_sles_hypre_t  *hypre_ctx = cs_sles_get_context(sles);
      assert(hypre_ctx);

      cs_sles_hypre_set_n_max_iter(hypre_ctx, cvgp.n_max_iter);
    }
    break;
#endif

  default:
    /* CS_PARAM_SOLVER_CLASS_MUMPS => Nothing to do */
    break;
  }
}

#if defined(HAVE_PETSC)
/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the command line option for PETSc
 *
 * \param[in] use_prefix     need a prefix
 * \param[in] prefix         optional prefix
 * \param[in] keyword        command keyword
 * \param[in] keyval         command value
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_setup_petsc_cmd(bool         use_prefix,
                              const char  *prefix,
                              const char  *keyword,
                              const char  *keyval)
{
  _petsc_cmd(use_prefix, prefix, keyword, keyval);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set a couple (preconditioner, solver) in PETSc
 *
 * \param[in]      label  label to identify this (part of) system
 * \param[in, out] slesp  pointer to a set of SLES parameters
 * \param[in, out] p_ksp  solver structure for PETSc
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_setup_petsc_ksp(const char       *label,
                              cs_param_sles_t  *slesp,
                              void             *p_ksp)
{
  KSP  ksp = p_ksp;
  assert(ksp != NULL);

  int len = strlen(label) + 1;
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
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set an AMG preconditioner in PETSc
 *
 * \param[in]      prefix  label to identify this (part of) system
 * \param[in]      slesp   pointer to a set of SLES parameters
 * \param[in, out] p_pc    preconditioner structure for PETsc
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_setup_petsc_pc_amg(const char       *prefix,
                                 cs_param_sles_t  *slesp,
                                 void             *p_pc)
{
  PC  pc = p_pc;
  assert(pc != NULL);
  bool  is_symm = _system_should_be_sym(slesp->solver);

  switch (slesp->amg_type) {

  case CS_PARAM_AMG_PETSC_GAMG_V:
  case CS_PARAM_AMG_PETSC_GAMG_W:
  case CS_PARAM_AMG_PETSC_PCMG:
    _petsc_pcgamg_hook(prefix, slesp, is_symm, pc);
    break;

  case CS_PARAM_AMG_HYPRE_BOOMER_V:
  case CS_PARAM_AMG_HYPRE_BOOMER_W:
    if (cs_param_sles_hypre_from_petsc())
      _petsc_pchypre_hook(prefix, slesp, is_symm, pc);

    else {

      cs_base_warn(__FILE__, __LINE__);
      cs_log_printf(CS_LOG_WARNINGS,
                    "%s: prefix=\"%s\": Switch to GAMG since BoomerAMG is not"
                    " available.\n",
                    __func__, prefix);
      _petsc_pcgamg_hook(prefix, slesp, is_symm, pc);

    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: prefix=\"%s\": Invalid AMG type for the PETSc library.",
              __func__, prefix);
    break;

  } /* End of switch on the AMG type */
}
#endif

/*----------------------------------------------------------------------------*/

END_C_DECLS