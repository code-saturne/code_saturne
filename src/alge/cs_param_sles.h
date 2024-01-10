#ifndef __CS_PARAM_SLES_H__
#define __CS_PARAM_SLES_H__

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_param_types.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
  \file cs_param_sles.h

  \brief Structure and routines handling the SLES settings stored inside a
         cs_param_sles_t structure

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

/* MUMPS settings */
/* ============== */

/*! \enum cs_param_sles_facto_type_t
 *  \brief type of factorization to consider when using the MUMPS solver to
 *  solve a linear system
 *
 * \var CS_PARAM_SLES_FACTO_LU
 * \brief LU factorization is the most generic factorization available with
 * MUMPS. It can handle general matrices (block and/or unsymmetric matrices)
 *
 * \var CS_PARAM_SLES_FACTO_LDLT_SYM
 * \brief This factorization is a Cholesky factorization (L.D.Lt) for general
 * symmetric matrices
 *
 * \var CS_PARAM_SLES_FACTO_LDLT_SPD
 * \brief This factorization is devoted to SPD matrices and corresponds to a
 * Cholesky factorization. This is more specific and thus more efficient than
 * \ref CS_PARAM_SLES_FACTO_LDLT_SYM
 */

typedef enum {

  CS_PARAM_SLES_FACTO_LU,
  CS_PARAM_SLES_FACTO_LDLT_SYM,
  CS_PARAM_SLES_FACTO_LDLT_SPD,

  CS_PARAM_SLES_N_FACTO_TYPES

} cs_param_sles_facto_type_t;

/*! \enum cs_param_sles_analysis_algo_t
 *  \brief Type of algorithm to consider when using the MUMPS solver to perform
 *  the analysis step (renumbering and graph manipulation). Please refer to the
 *  MUMPS user guide for more details about the following algorithms. AMD, QAMD
 *  and PORD are available with MUMPS without any prerequesite.
 *
 * \var CS_PARAM_SLES_ANALYSIS_AMD
 * AMD is a sequential algorithm which is well-suited for 2D problem (for 3D
 * problems it induces a higher memeory consumption).
 *
 * \var CS_PARAM_SLES_ANALYSIS_QAMD
 * QAMD is a sequential algorithm which is well-suited for 2D problem (for 3D
 * problems it induces a higher memeory consumption).
 *
 * \var CS_PARAM_SLES_ANALYSIS_PORD
 * PORD is a sequential algorithm which is a good trade-off when MUMPS is
 * installed with no prerequisite such as METIS or Scotch.
 *
 * \var CS_PARAM_SLES_ANALYSIS_SCOTCH
 * SCOTCH is a sequential algorithm which delivers the very good performance
 * with 3D meshes, generally, better than PORD and not as good as METIS
 *
 * \var CS_PARAM_SLES_ANALYSIS_PTSCOTCH
 * PTSCOTCH is a parallel version of the sequential SCOTCH algorithm
 *
 * \var CS_PARAM_SLES_ANALYSIS_METIS
 * METIS is a sequential algorithm which delivers the best performance in case
 * of 2D meshes.
 *
 * \var CS_PARAM_SLES_ANALYSIS_PARMETIS
 * PARMETIS is a parallel version of the sequential METIS algorithm
 *
 * \var CS_PARAM_SLES_ANALYSIS_AUTO
 * MUMPS decides what is the best choice among available algorithms. This is
 * the default choice.
 */

typedef enum {

  CS_PARAM_SLES_ANALYSIS_AMD,
  CS_PARAM_SLES_ANALYSIS_QAMD,
  CS_PARAM_SLES_ANALYSIS_PORD,
  CS_PARAM_SLES_ANALYSIS_SCOTCH,
  CS_PARAM_SLES_ANALYSIS_PTSCOTCH,
  CS_PARAM_SLES_ANALYSIS_METIS,
  CS_PARAM_SLES_ANALYSIS_PARMETIS,

  CS_PARAM_SLES_ANALYSIS_AUTO,

  CS_PARAM_SLES_N_ANALYSIS_ALGOS

} cs_param_sles_analysis_algo_t;


/*! \enum cs_param_sles_memory_usage_t
 *  \brief Strategy for the memory usage inside MUMPS
 *
 * \var CS_PARAM_SLES_MEMORY_CONSTRAINED
 * Strategy aiming at limiting the memory usage
 *
 * \var CS_PARAM_SLES_MEMORY_AUTO
 * Strategy relying on the default settings
 *
 * \var CS_PARAM_SLES_MEMORY_CPU_DRIVEN
 * Strategy aiming at the best CPU time
 */

typedef enum {

  CS_PARAM_SLES_MEMORY_CONSTRAINED,
  CS_PARAM_SLES_MEMORY_AUTO,
  CS_PARAM_SLES_MEMORY_CPU_DRIVEN,

  CS_PARAM_SLES_N_MEMORY_USAGES

} cs_param_sles_memory_usage_t;

/*! \struct cs_param_sles_mumps_t
 *  \brief Set of parameters to specify additional options to MUMPS
 *  For more advanced settings, one has to use the \ref cs_user_sles_mumps_hook
 *  function. Please also refer to the MUMPS user guide for more details.
 */

typedef struct {

  /* \var analysis_algo
   * Choice of the algorithm used to perform the analysis step
   *
   * \var facto_type
   * Type of factorization to consider. This choice depends on the type of
   * matrix to handle
   *
   * \var mem_usage
   * Type of strategy to consider for the memory usage.
   */

  cs_param_sles_analysis_algo_t   analysis_algo;
  cs_param_sles_facto_type_t      facto_type;
  cs_param_sles_memory_usage_t    mem_usage;

  bool    is_single;        /*!< Single precision is used, otherwise double */

  bool    advanced_optim;   /*!< Activate advanced optimizations (very useful
                                 when openMP is used) */

  double  blr_threshold;    /*!< Dropping parameter in the BLR compression. The
                                 value is directly related to the accuracy of
                                 the compression (0: not used). A positive
                                 value implies a usage of the predefined
                                 algorithm. A negative value implies the usage
                                 of an alternative algorithm. */

  double  mem_coef;         /*!< Percentage of increase of the
                                 automatically-defined memory workspace. Really
                                 useful for 2D cases. (Not used if < 0) */

  int     block_analysis;   /*!< Analysis is performed by block. Value of the
                                 block size. Not used if < 1 */

  int     ir_steps;         /*!< Number of steps for the Iterative Refinement */

} cs_param_sles_mumps_t;

/* BoomerAMG settings */
/* ================== */

/*! \enum cs_param_sles_boomeramg_coarsen_algo_t
 *  \brief Type of algorithm used in boomerAMG to coarsen a level. Only a
 *  selection of algorithms is available here. Values are those given in HYPRE
 */

typedef enum {

  CS_PARAM_SLES_BOOMERAMG_COARSEN_FALGOUT = 6,
  CS_PARAM_SLES_BOOMERAMG_COARSEN_PMIS = 8,
  CS_PARAM_SLES_BOOMERAMG_COARSEN_HMIS = 10, /* (default) */
  CS_PARAM_SLES_BOOMERAMG_COARSEN_CGC = 21,
  CS_PARAM_SLES_BOOMERAMG_COARSEN_CGC_E = 22,

  CS_PARAM_SLES_BOOMERAMG_N_COARSEN_ALGOS

} cs_param_sles_boomeramg_coarsen_algo_t;

/*! \enum cs_param_sles_boomeramg_interp_type_t
 *  \brief Type of algorithm used in boomerAMG to coarsen a level. Only a
 *  selection of algorithms is available here. Values are those given in HYPRE
 */

typedef enum {

  CS_PARAM_SLES_BOOMERAMG_INTERP_HYPERBOLIC = 2,
  CS_PARAM_SLES_BOOMERAMG_INTERP_EXT_PLUS_I_CC = 6, /* (default) Also for GPU */
  CS_PARAM_SLES_BOOMERAMG_INTERP_EXT_PLUS_I = 7,
  CS_PARAM_SLES_BOOMERAMG_INTERP_FF1 = 13,
  CS_PARAM_SLES_BOOMERAMG_INTERP_EXTENDED = 14,     /* Also for GPU */
  CS_PARAM_SLES_BOOMERAMG_INTERP_EXT_PLUS_I_MATRIX = 17,
  CS_PARAM_SLES_BOOMERAMG_INTERP_EXT_PLUS_E_MATRIX = 18,

  CS_PARAM_SLES_BOOMERAMG_N_INTERP_ALGOS

} cs_param_sles_boomeramg_interp_algo_t;

/*! \enum cs_param_sles_boomeramg_interp_type_t
 *  \brief Type of algorithm used in boomerAMG to coarsen a level. Only a
 *  selection of algorithms is available here. Values are those used in HYPRE.
 */

typedef enum {

  CS_PARAM_SLES_BOOMERAMG_JACOBI = 0,
  CS_PARAM_SLES_BOOMERAMG_FORWARD_GS = 3,
  CS_PARAM_SLES_BOOMERAMG_BACKWARD_GS = 4,
  CS_PARAM_SLES_BOOMERAMG_HYBRID_SSOR = 6,
  CS_PARAM_SLES_BOOMERAMG_L1_SGS = 8,
  CS_PARAM_SLES_BOOMERAMG_GAUSS_ELIM = 9,      /* for the coarsest level only */
  CS_PARAM_SLES_BOOMERAMG_BACKWARD_L1_GS = 13, /* (default) */
  CS_PARAM_SLES_BOOMERAMG_FORWARD_L1_GS = 14,  /* (default) */
  CS_PARAM_SLES_BOOMERAMG_CG = 15,
  CS_PARAM_SLES_BOOMERAMG_CHEBYSHEV = 16,
  CS_PARAM_SLES_BOOMERAMG_FCF_JACOBI = 17,
  CS_PARAM_SLES_BOOMERAMG_L1_JACOBI = 18,

  CS_PARAM_SLES_BOOMERAMG_N_SMOOTHERS

} cs_param_sles_boomeramg_smoother_t;

/*! \struct cs_param_sles_boomeramg_boomerang_t
 *  \brief Set of the most influencial parameters to setup the algebraic
 *   multigrid BoomerAMG belonging to the HYPRE library. These parameters are
 *   used to define this AMG directly in HYPRE or through the PETSc library
 *   according to the settings and the installed dependencies. Please refer to
 *  the HYPRE documentation for more details.
 */

typedef struct {

  /* Read the function \ref _petsc_pchypre_hook or \ref _hypre_boomeramg_hook
     for more details and read the HYPRE user guide */

  double                                   strong_threshold;
  cs_param_sles_boomeramg_coarsen_algo_t   coarsen_algo;
  cs_param_sles_boomeramg_interp_algo_t    interp_algo;
  int                                      p_max;
  int                                      n_agg_levels;
  int                                      n_agg_paths;

  cs_param_sles_boomeramg_smoother_t       down_smoother;
  cs_param_sles_boomeramg_smoother_t       up_smoother;
  cs_param_sles_boomeramg_smoother_t       coarse_solver;

  int                                      n_down_iter;
  int                                      n_up_iter;

} cs_param_sles_boomeramg_t;

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
   *  systems or more complex systems */

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

/*!
 * \struct cs_param_sles_saddle_t
 * \brief Structure storing all metadata related to the resolution of a
 *        saddle-point linear system
 */

typedef struct {

  /*!
    * \var verbosity
    * verbosity (level of information displayed)
    */

  int                         verbosity;

  /*!
   * @name Main parameter settings
   *
   * Set of parameters to drive the resolution of a saddle-point system
   */
  /*! @{ */

  /*!
   * \var solver
   *  Type of solver to solve the saddle-point system
   *  If solver is set to CS_PARAM_SADDLE_N_SOLVERS, then there is no need to
   *  solve a saddle-point system.
   */

  cs_param_saddle_solver_t    solver;

  /*! \var precond
   *  Type of preconditioner for the saddle-point system which is viewed as a
   *  2x2 block matrix.
   */

  cs_param_saddle_precond_t   precond;

  /*! \var cvg_param
   *  Structure storing the parameters to know if an iterative process has to
   *  stop (convergence or divergence). These criteria are related the
   *  iterative algorithm used to solve the saddle-point system. This is the
   *  case for instance with an Uzawa or GKB algorithm.
   */

  cs_param_sles_cvg_t         cvg_param;

  /*! @} */

  /*!
   * @name Schur complement approximation
   *
   * Set of parameters to drive the resolution of the pressure-related
   * block. This is often a Schur complement approximation to B.A^-1.Bt
   */
  /*! @{ */

  /*!
   * \var schur_approximation
   *
   * \var schur_approximation
   * Choice of the way of preconditioning the schur approximation
   */

  cs_param_schur_approx_t     schur_approximation;

  /*! \var schur_sles_param
   * Set of paremeters used to solve the Schur complement if needed. This
   * depends on the type of Schur approximation which has been cosen.
   */

  cs_param_sles_t            *schur_sles_param;

  /*! @} */

} cs_param_sles_saddle_t;

/*============================================================================
 * Global variables
 *============================================================================*/

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
cs_param_sles_saddle_create(void);

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
                                cs_param_sles_saddle_t    *saddlep);

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
                          cs_param_sles_saddle_t         *dest);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free the structure storing the parameter settings for a saddle-point
 *        system
 *
 * \param[in, out] p_saddlep    double pointer to the structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_saddle_free(cs_param_sles_saddle_t    **p_saddlep);

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
 * \brief Set the main members of a cs_param_sles_mumps_t structure. This
 *        structure is allocated if needed. Other members are kept to their
 *        values.
 *
 * \param[in, out] slesp           pointer to a cs_param_sles_t structure
 * \param[in]      coarsen_algo    type of algoritmh for the coarsening
 * \param[in]      down_smoother   type of smoother for the down cycle
 * \param[in]      up_smoother     type of smoother for th up cycle
 * \param[in]      coarse_solver   solver at the coarsest level
 * \param[in]      n_down_iter     number of smoothing steps for the down cycle
 * \param[in]      n_up_iter       number of smoothing steps for the up cycle
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_boomeramg(cs_param_sles_t                          *slesp,
                        cs_param_sles_boomeramg_coarsen_algo_t    coarsen_algo,
                        cs_param_sles_boomeramg_smoother_t        down_smoother,
                        cs_param_sles_boomeramg_smoother_t        up_smoother,
                        cs_param_sles_boomeramg_smoother_t        coarse_solver,
                        int                                       n_down_iter,
                        int                                       n_up_iter);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the members of a cs_param_sles_boomeramg_t structure used in
 *        advanced settings. This structure is allocated if needed. Other
 *        members are kept to their values. Please refer to the HYPRE user
 *        guide for more details about the following advanced options.
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
cs_param_sles_boomeramg_advanced(cs_param_sles_t                  *slesp,
                                 double                            strong_thr,
                             cs_param_sles_boomeramg_interp_algo_t interp_algo,
                                 int                               p_max,
                                 int                               n_agg_lv,
                                 int                               n_agg_paths);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the main members of a cs_param_sles_mumps_t structure. This
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
                    cs_param_sles_facto_type_t    facto_type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the members of a cs_param_sles_mumps_t structure used in advanced
 *        settings. This structure is allocated if needed. Other members are
 *        kept to their values. Please refer to the MUMPS user guide for more
 *        details about the following advanced options.
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
                             bool                            advanced_optim);

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
 * \param[in]      use_field_id  if false use system name to define a SLES
 * \param[in, out] slesp         pointer to a cs_param_sles_t structure
 *
 * \return an error code (-1 if a problem is encountered, 0 otherwise)
 */
/*----------------------------------------------------------------------------*/

int
cs_param_sles_set(bool                 use_field_id,
                  cs_param_sles_t     *slesp);

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
 * \brief Retrieve the related solver class from the amg type
 *
 * \param[in]  amg_type    type of AMG to consider
 *
 * \return the related solver class or CS_PARAM_SLES_CLASS_CS
 */
/*----------------------------------------------------------------------------*/

cs_param_sles_class_t
cs_param_sles_get_class_from_amg(cs_param_amg_type_t   amg_type);

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

/*============================================================================
 * Static inline public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Find if a MUMPS-related solver is set or not
 *
 * \param[in] solver   type of solver
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

static inline bool
cs_param_sles_is_mumps_set(cs_param_itsol_type_t  solver)
{
  if (solver == CS_PARAM_ITSOL_MUMPS)
    return true;
  else
    return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a default context storing a set of BoomerAMG parameters
 *
 * \param[in] slesp     pointer to a cs_sles_param_t structure to update
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_param_sles_set_default_boomeramg_context(cs_param_sles_t  *slesp)
{
  assert(slesp != NULL);
  BFT_FREE(slesp->context_param); /* Avoid an issue if the context was
                                     different */

  cs_param_sles_boomeramg(slesp,
                          CS_PARAM_SLES_BOOMERAMG_COARSEN_HMIS,
                          CS_PARAM_SLES_BOOMERAMG_HYBRID_SSOR,
                          CS_PARAM_SLES_BOOMERAMG_HYBRID_SSOR,
                          CS_PARAM_SLES_BOOMERAMG_GAUSS_ELIM,
                          1,  /* n_down_iter */
                          1); /* n_up_iter */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PARAM_SLES_H__ */
