#ifndef __CS_PARAM_SLES_H__
#define __CS_PARAM_SLES_H__

/*============================================================================
 * Routines to handle the SLES (Sparse Linear Equation Solver) settings
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_param_amg.h"
#include "cs_param_mumps.h"
#include "cs_param_types.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

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

/*! \struct cs_param_sles_cvg_t
 *  \brief Set of parameters to check the convergence (or the divergence) of an
 *         iterative process (tolerances or max. number of iterations)
 */

typedef struct {

/*!
 * \var atol
 * Absolute tolerance under which the iterative process is stopped
 *
 * \var rtol
 * Relative tolerance under which the iterative process is stopped
 *
 * \var dtol
 * Tolerance above which the iterative process is stated as "diverged".
 * Not used if < 0
 *
 * \var n_max_iter
 * Maximal number of iterations before stopping the iterative process
 */

  double               atol;
  double               rtol;
  double               dtol;
  int                  n_max_iter;

} cs_param_sles_cvg_t;

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
  bool    setup_done;  /*!< SLES setup step has been done */

  cs_param_sles_class_t      solver_class; /*!< class of SLES to consider  */
  cs_param_precond_type_t    precond;      /*!< type of preconditioner */
  cs_param_itsol_type_t      solver;       /*!< type of solver */
  bool                       flexible;     /*!< need a flexible variant ? */
  int                        restart;      /*!< max. iter. before restarting */
  cs_param_amg_type_t        amg_type;     /*!< type of AMG algorithm */

  /*! \var pcd_block_type
   *  type of block preconditioner to use (only meaningful for vector-valued
   *  systems or more complex systems) */

  cs_param_precond_block_t    pcd_block_type;

  /*! \var resnorm_type
   *  normalized or not the norm of the residual used for the stopping criterion
   *  See \ref CS_EQKEY_ITSOL_RESNORM_TYPE for more details. */

  cs_param_resnorm_type_t     resnorm_type;

  /*! \var cvg_param
   *  Structure storing the parameters to know if an iterative process has to
   *  stop (convergence or divergence).
   */

  cs_param_sles_cvg_t         cvg_param;

  /*! \var context_param
   *  Pointer to a structure cast on-the-fly storing specific parameters
   *  related to the given solver (MUMPS for instance)
   */

  void                       *context_param;

} cs_param_sles_t;

/*============================================================================
 * Global variables
 *============================================================================*/

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
cs_param_sles_create(int           field_id,
                     const char   *system_name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free a \ref cs_param_sles_t structure
 *
 * \param[in, out] slesp    pointer to a \cs_param_sles_t structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_free(cs_param_sles_t   **p_slesp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log information related to the linear settings stored in the
 *        structure
 *
 * \param[in] slesp    pointer to a \ref cs_param_sles_log
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_log(cs_param_sles_t   *slesp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Copy a cs_param_sles_t structure from src to dst
 *
 * \param[in]      src    reference cs_param_sles_t structure to copy
 * \param[in, out] dst    copy of the reference at exit
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_copy_from(const cs_param_sles_t   *src,
                        cs_param_sles_t         *dst);

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
                  cs_param_sles_t     *slesp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate and initialize a new context structure for the boomerAMG
 *        settings.
 *
 * \param[in, out] slesp        pointer to a cs_param_sles_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_boomeramg_reset(cs_param_sles_t  *slesp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the main members of a cs_param_amg_boomer_t structure. This
 *        structure is allocated and initialized and then one sets the main
 *        given parameters. Please refer to the HYPRE user guide for more
 *        details about the following options.
 *
 * \param[in, out] slesp           pointer to a cs_param_sles_t structure
 * \param[in]      n_down_iter     number of smoothing steps for the down cycle
 * \param[in]      down_smoother   type of smoother for the down cycle
 * \param[in]      n_up_iter       number of smoothing steps for the up cycle
 * \param[in]      up_smoother     type of smoother for th up cycle
 * \param[in]      coarse_solver   solver at the coarsest level
 * \param[in]      coarsen_algo    type of algoritmh for the coarsening
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
/*!
 * \brief Set the members of a cs_param_amg_boomer_t structure used in
 *        advanced settings. This structure is allocated if needed. Other
 *        members are kept to their values. Please refer to the HYPRE user
 *        guide for more details about the following options.
 *
 * \param[in, out] slesp            pointer to a cs_param_sles_t structure
 * \param[in]      strong_thr       value of the strong threshold (coarsening)
 * \param[in]      interp_algo      algorithm used for the interpolation
 * \param[in]      p_max            max number of elements per row (interp)
 * \param[in]      n_agg_lv         aggressive coarsening (number of levels)
 * \param[in]      n_agg_paths      aggressive coarsening (number of paths)
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_boomeramg_advanced(cs_param_sles_t                   *slesp,
                                 double                             strong_thr,
                                 cs_param_amg_boomer_interp_algo_t  interp_algo,
                                 int                                p_max,
                                 int                                n_agg_lv,
                                 int                                n_agg_paths);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the main members of a cs_param_mumps_t structure. This structure
 *        is allocated and initialized with default settings if needed. If the
 *        structure exists already, then advanced members are kept to their
 *        values.
 *
 * \param[in, out] slesp         pointer to a cs_param_sles_t structure
 * \param[in]      is_single     single-precision or double-precision
 * \param[in]      facto_type    type of factorization to consider
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_mumps(cs_param_sles_t             *slesp,
                    bool                         is_single,
                    cs_param_mumps_facto_type_t  facto_type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the members related to an advanced settings of a cs_param_mumps_t
 *        structure. This structure is allocated and initialized if
 *        needed. Please refer to the MUMPS user guide for more details about
 *        the following advanced options.
 *
 * \param[in, out] slesp            pointer to a cs_param_sles_t structure
 * \param[in]      analysis_algo    algorithm used for the analysis step
 * \param[in]      block_analysis   > 1: fixed block size; otherwise do nothing
 * \param[in]      mem_coef         percentage increase in the memory workspace
 * \param[in]      blr_threshold    Accuracy in BLR compression (0: not used)
 * \param[in]      ir_steps         0: No, otherwise the number of iterations
 * \param[in]      mem_usage        strategy to adopt for the memory usage
 * \param[in]      advanced_optim   activate advanced optimization (MPI/openMP)
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_mumps_advanced(cs_param_sles_t                *slesp,
                             cs_param_mumps_analysis_algo_t  analysis_algo,
                             int                             block_analysis,
                             double                          mem_coef,
                             double                          blr_threshold,
                             int                             ir_steps,
                             cs_param_mumps_memory_usage_t   mem_usage,
                             bool                            advanced_optim);

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
                                  const cs_param_sles_t   *slesp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check the availability of Hypre solvers from the PETSc library
 *
 * \return return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_param_sles_hypre_from_petsc(void);

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
cs_param_sles_check_class(cs_param_sles_class_t   wanted_class);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check if the setting related to the AMG is consistent with the
 *        solver class. If an issue is detected, try to solve it whith the
 *        nearest option.
 *
 * \param[in, out] slesp    pointer to a cs_param_sles_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_check_amg(cs_param_sles_t   *slesp);

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
                        const char   *keyval);
#endif

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PARAM_SLES_H__ */
