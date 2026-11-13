#ifndef CS_PARAM_SLES_H
#define CS_PARAM_SLES_H

/*============================================================================
 * Routines to handle the SLES (Sparse Linear Equation Solver) settings
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "alge/cs_param_amg.h"
#include "alge/cs_param_hpddm.h"
#include "alge/cs_param_mumps.h"
#include "base/cs_param_types.h"

/*----------------------------------------------------------------------------*/

/*!
  \file cs_param_sles.h

  \brief Structure and routines handling the SLES ((Sparse Linear Equation
         Solver) settings stored inside a cs_param_sles_t structure

*/

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*!
 * \struct cs_param_sles_t
 * \brief Structure storing all metadata related to the resolution of a linear
 *        system with an iterative solver.
 */

typedef struct {

  char   *name;        /*!< System name (equation name if this is automatic) */
  int     field_id;    /*!< Field id related to a SLES. By default, this is set
                            to -1 */
  int     verbosity;   /*!< SLES verbosity */

  cs_param_solver_class_t    solver_class;  /*!< class of SLES to consider  */
  cs_param_precond_type_t    precond;       /*!< type of preconditioner */
  cs_param_solver_type_t     solver;        /*!< type of solver */
  bool                       need_flexible; /*!< need a flexible variant ? */
  int                        restart;       /*!< max. iter. before restarting */
  cs_param_amg_type_t        amg_type;      /*!< type of AMG algorithm */

  /*! \var precond_block_type
   *  type of block preconditioner to use (only meaningful for vector-valued
   *  systems or more complex systems) */

  cs_param_precond_block_t   precond_block_type;

  /*! \var resnorm_type
   *  normalized or not the norm of the residual used for the stopping criterion
   *  See \ref CS_EQKEY_SOLVER_RESNORM_TYPE for more details. */

  cs_param_resnorm_type_t    resnorm_type;

  /*! \var allow_no_op
   *  allow one to not call the solve function in some specific situations
   *  Please refer to \c _needs_solving for more details
   */

  bool                       allow_no_op;

  /*! \var mat_is_sym
   *  allow one to know if the matric is symmetric (independent of the way the
   *  storage is done).
   * Default value is false
   */

  bool                       mat_is_sym;

  /*! \var cvg_param
   *  Structure storing the parameters to know if an iterative process has to
   *  stop (convergence or divergence).
   */

  cs_param_convergence_t     cvg_param;

  /*! \var context_param
   *  Pointer to a structure cast on-the-fly storing specific parameters
   *  related to the given solver (MUMPS for instance)
   */

  void                      *context_param;

} cs_param_sles_t;

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * Create a cs_param_sles_t structure and assign a default settings
 */
/*----------------------------------------------------------------------------*/

cs_param_sles_t *
cs_param_sles_create(int           field_id,
                     const char   *system_name);

/*----------------------------------------------------------------------------*/
/*
 * Free a cs_param_sles_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_free(cs_param_sles_t  **p_slesp);

/*----------------------------------------------------------------------------*/
/*
 * Log information related to the linear settings stored in the structure
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_log(cs_param_sles_t   *slesp);

/*----------------------------------------------------------------------------*/
/*
 * Copy a cs_param_sles_t structure from src to dst.
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_copy_from(const cs_param_sles_t  *src,
                        cs_param_sles_t        *dst);

/*----------------------------------------------------------------------------*/
/*
 * Set the solver associated to this SLES from its keyval.
 */
/*----------------------------------------------------------------------------*/

int
cs_param_sles_set_solver(const char       *keyval,
                         cs_param_sles_t  *slesp);

/*----------------------------------------------------------------------------*/
/*
 * Set the preconditioner associated to this SLES from its keyval.
 */
/*----------------------------------------------------------------------------*/

int
cs_param_sles_set_precond(const char       *keyval,
                          cs_param_sles_t  *slesp);

/*----------------------------------------------------------------------------*/
/*
 * Set the class of solvers associated to this SLES from its keyval.
 * Common choices are "petsc", "hypre" or "mumps".
 */
/*----------------------------------------------------------------------------*/

int
cs_param_sles_set_solver_class(const char       *keyval,
                               cs_param_sles_t  *slesp);

/*----------------------------------------------------------------------------*/
/*
 * Set the type of algebraic multigrid (AMG) associated to this SLES from its
 * keyval.
 */
/*----------------------------------------------------------------------------*/

int
cs_param_sles_set_amg_type(const char       *keyval,
                           cs_param_sles_t  *slesp);

/*----------------------------------------------------------------------------*/
/*
 * Set the type of block preconditioner associated to this SLES from its
 * keyval.
 */
/*----------------------------------------------------------------------------*/

int
cs_param_sles_set_precond_block_type(const char       *keyval,
                                     cs_param_sles_t  *slesp);

/*----------------------------------------------------------------------------*/
/*
 * Set the convergence criteria for the given SLES parameters. One can use the
 * parameter value CS_CDO_KEEP_DEFAULT to let unchange the settings.
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_set_cvg_param(cs_param_sles_t  *slesp,
                            double            rtol,
                            double            atol,
                            double            dtol,
                            int               max_iter);

/*----------------------------------------------------------------------------*/
/*
 * Allocate and initialize a new context structure for the in-house AMG
 * settings.
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_amg_inhouse_reset(cs_param_sles_t  *slesp,
                                bool              used_as_solver,
                                bool              used_as_k_cycle);

/*----------------------------------------------------------------------------*/
/*
 * Set the main members of a cs_param_amg_inhouse_t structure.
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_amg_inhouse(cs_param_sles_t                *slesp,
                          int                             n_down_iter,
                          cs_param_amg_inhouse_solver_t   down_smoother,
                          int                             down_poly_deg,
                          int                             n_up_iter,
                          cs_param_amg_inhouse_solver_t   up_smoother,
                          int                             up_poly_deg,
                          cs_param_amg_inhouse_solver_t   coarse_solver,
                          int                             coarse_poly_deg,
                          cs_param_amg_inhouse_coarsen_t  coarsen_algo,
                          int                             aggreg_limit);

/*----------------------------------------------------------------------------*/
/*
 * Set the members of a cs_param_amg_inhouse_t structure used in advanced
 * settings.
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_amg_inhouse_advanced(cs_param_sles_t  *slesp,
                                   int               max_levels,
                                   cs_gnum_t         min_n_g_rows,
                                   double            p0p1_relax,
                                   int               coarse_max_iter,
                                   double            coarse_rtol_mult,
                                   int               fine_lvl_agg,
                                   int               fine_coarsening);

/*----------------------------------------------------------------------------*/
/*
 * Allocate and initialize a new context structure for the boomerAMG settings.
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_boomeramg_reset(cs_param_sles_t  *slesp);

/*----------------------------------------------------------------------------*/
/*
 * Set the main members of a cs_param_amg_boomer_t structure.
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_boomeramg(cs_param_sles_t                    *slesp,
                        int                                 n_down_iter,
                        cs_param_amg_boomer_smoother_t      down_smoother,
                        int                                 n_up_iter,
                        cs_param_amg_boomer_smoother_t      up_smoother,
                        cs_param_amg_boomer_smoother_t      coarse_solver,
                        cs_param_amg_boomer_coarsen_algo_t  coarsen_algo);

/*----------------------------------------------------------------------------*/
/*
 * Set the members of a cs_param_amg_boomer_t structure used in advanced
 * settings.
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_boomeramg_advanced(cs_param_sles_t                  *slesp,
                                 double                            strong_thr,
                                 cs_param_amg_boomer_interp_algo_t interp_algo,
                                 int                               p_max,
                                 int                               n_agg_lv,
                                 int                               n_agg_paths);

/*----------------------------------------------------------------------------*/
/*
 * Allocate and initialize a new context structure for GAMG in PETSc.
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_gamg_reset(cs_param_sles_t  *slesp);

/*----------------------------------------------------------------------------*/
/*
 * Set the main members of a cs_param_amg_gamg_t structure.
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_gamg(cs_param_sles_t                  *slesp,
                   int                               n_down_iter,
                   cs_param_amg_gamg_smoother_t      down_smoother,
                   int                               n_up_iter,
                   cs_param_amg_gamg_smoother_t      up_smoother,
                   cs_param_amg_gamg_coarse_solver_t coarse_solver);

/*----------------------------------------------------------------------------*/
/*
 * Set the members of a cs_param_amg_gamg_t structure used in advanced settings.
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_gamg_advanced(cs_param_sles_t *slesp,
                            double           threshold,
                            int              n_agg_lv,
                            bool             use_sq_grph,
                            int              n_smooth_agg);

/*----------------------------------------------------------------------------*/
/*
 * Allocate and initialize a new context structure for HMG in PETSc.
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_hmg_reset(cs_param_sles_t  *slesp);

/*----------------------------------------------------------------------------*/
/*
 * Set the main members of a cs_param_amg_hmg_t structure.
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_hmg(cs_param_sles_t                   *slesp,
                  int                                n_down_iter,
                  cs_param_amg_gamg_smoother_t       down_smoother,
                  int                                n_up_iter,
                  cs_param_amg_gamg_smoother_t       up_smoother,
                  cs_param_amg_gamg_coarse_solver_t  coarse_solver);

/*----------------------------------------------------------------------------*/
/*
 * Set the members of a cs_param_amg_gamg_t structure used in advanced settings.
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_hmg_advanced(cs_param_sles_t *slesp,
                           bool             use_boomer_coarsening,
                           bool             reuse_interpolation,
                           bool             subspace_coarsening);

/*----------------------------------------------------------------------------*/
/*
 * Allocate and initialize a new context structure for MUMPS settings.
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_mumps_reset(cs_param_sles_t  *slesp);

/*----------------------------------------------------------------------------*/
/*
 * Set the main members of a cs_param_mumps_t structure.
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_mumps(cs_param_sles_t             *slesp,
                    bool                         is_single,
                    cs_param_mumps_facto_type_t  facto_type);

/*----------------------------------------------------------------------------*/
/*
 * Set the members related to an advanced settings of a cs_param_mumps_t structure.
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_mumps_advanced(cs_param_sles_t                *slesp,
                             cs_param_mumps_analysis_algo_t  analysis_algo,
                             int                             block_analysis,
                             bool                            keep_ordering,
                             double                          mem_coef,
                             double                          blr_threshold,
                             int                             ir_steps,
                             cs_param_mumps_memory_usage_t   mem_usage,
                             bool                            advanced_optim);

/*----------------------------------------------------------------------------*/
/*
 * Set keep_data option for MUMPS solver in SLES parameters.
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_mumps_set_keep_data(cs_param_sles_t  *slesp,
                                  bool              keep_data);

/*----------------------------------------------------------------------------*/
/*
 * Query keep_data option for MUMPS solver in SLES parameters.
 */
/*----------------------------------------------------------------------------*/

bool
cs_param_sles_mumps_get_keep_data(const cs_param_sles_t  *slesp);

/*----------------------------------------------------------------------------*/
/*
 * Allocate and initialize a new context structure for HPDDM settings.
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_hpddm_reset(cs_param_sles_t *slesp);

/*----------------------------------------------------------------------------*/
/*
 * Set the members related to an settings of a cs_param_hpddm_t structure.
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_hpddm(cs_param_sles_t *slesp,
                    const bool       use_neumann,
                    const int        nb_eigenvector,
                    const int        harmonic_overlap,
                    const double     relative_threshold,
                    const int        p = -1);

/*----------------------------------------------------------------------------*/
/*
 * Set the members related to an advanced settings of a cs_param_hpddm_t
 * structure.
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_hpddm_advanced(cs_param_sles_t *slesp,
                             const int        min_iter,
                             const int        max_iter,
                             const int        min_harmonic_overlap,
                             const int        max_harmonic_overlap,
                             const int        min_nb_eigenvector,
                             const int        max_nb_eigenvector);

/*----------------------------------------------------------------------------*/
/*
 * Check the availability of Hypre solvers from the PETSc library
 */
/*----------------------------------------------------------------------------*/

bool
cs_param_sles_hypre_from_petsc(void);

/*----------------------------------------------------------------------------*/
/*
 * Check the availability of a solver library and return the requested one if
 * this is possible or an alternative or CS_PARAM_N_SOLVER_CLASSES if no
 * alternative is available.
 */
/*----------------------------------------------------------------------------*/

cs_param_solver_class_t
cs_param_sles_check_class(cs_param_solver_class_t  wanted_class);

/*----------------------------------------------------------------------------*/

#endif /* CS_PARAM_SLES_H */
