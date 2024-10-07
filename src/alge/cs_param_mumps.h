#ifndef __CS_PARAM_MUMPS_H__
#define __CS_PARAM_MUMPS_H__

/*============================================================================
 * Routines and structure to handle the MUMPS settings
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

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
  \file cs_param_mumps.h

  \brief Routines and structure to handle the MUMPS setup. The structure is
         used as a context structure of a \ref cs_param_sles_t structure
*/

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! \enum cs_param_mumps_facto_type_t
 *  \brief type of factorization to consider when using the MUMPS solver to
 *  solve a linear system
 *
 * \var CS_PARAM_MUMPS_FACTO_LU
 * \brief LU factorization is the most generic factorization available with
 * MUMPS. It can handle general matrices (block and/or unsymmetric matrices)
 *
 * \var CS_PARAM_MUMPS_FACTO_LDLT_SYM
 * \brief This factorization is a Cholesky factorization (L.D.Lt) for general
 * symmetric matrices
 *
 * \var CS_PARAM_MUMPS_FACTO_LDLT_SPD
 * \brief This factorization is devoted to SPD matrices and corresponds to a
 * Cholesky factorization. This is more specific and thus more efficient than
 * \ref CS_PARAM_MUMPS_FACTO_LDLT_SYM
 */

typedef enum {

  CS_PARAM_MUMPS_FACTO_LU,
  CS_PARAM_MUMPS_FACTO_LDLT_SYM,
  CS_PARAM_MUMPS_FACTO_LDLT_SPD,

  CS_PARAM_MUMPS_N_FACTO_TYPES

} cs_param_mumps_facto_type_t;


/*! \enum cs_param_mumps_analysis_algo_t
 *  \brief Type of algorithm to consider when using the MUMPS solver to perform
 *  the analysis step (renumbering and graph manipulation). Please refer to the
 *  MUMPS user guide for more details about the following algorithms. AMD, QAMD
 *  and PORD are available with MUMPS without any prerequesite.
 *
 * \var CS_PARAM_MUMPS_ANALYSIS_AMD
 * AMD is a sequential algorithm which is well-suited for 2D problem (for 3D
 * problems it induces a higher memeory consumption).
 *
 * \var CS_PARAM_MUMPS_ANALYSIS_QAMD
 * QAMD is a sequential algorithm which is well-suited for 2D problem (for 3D
 * problems it induces a higher memeory consumption).
 *
 * \var CS_PARAM_MUMPS_ANALYSIS_PORD
 * PORD is a sequential algorithm which is a good trade-off when MUMPS is
 * installed with no prerequisite such as METIS or Scotch.
 *
 * \var CS_PARAM_MUMPS_ANALYSIS_SCOTCH
 * SCOTCH is a sequential algorithm which delivers the very good performance
 * with 3D meshes, generally, better than PORD and not as good as METIS
 *
 * \var CS_PARAM_MUMPS_ANALYSIS_PTSCOTCH
 * PTSCOTCH is a parallel version of the sequential SCOTCH algorithm
 *
 * \var CS_PARAM_MUMPS_ANALYSIS_METIS
 * METIS is a sequential algorithm which delivers the best performance in case
 * of 2D meshes.
 *
 * \var CS_PARAM_MUMPS_ANALYSIS_PARMETIS
 * PARMETIS is a parallel version of the sequential METIS algorithm
 *
 * \var CS_PARAM_MUMPS_ANALYSIS_AUTO
 * MUMPS decides what is the best choice among available algorithms. This is
 * the default choice.
 */

typedef enum {

  CS_PARAM_MUMPS_ANALYSIS_AMD,
  CS_PARAM_MUMPS_ANALYSIS_QAMD,
  CS_PARAM_MUMPS_ANALYSIS_PORD,
  CS_PARAM_MUMPS_ANALYSIS_SCOTCH,
  CS_PARAM_MUMPS_ANALYSIS_PTSCOTCH,
  CS_PARAM_MUMPS_ANALYSIS_METIS,
  CS_PARAM_MUMPS_ANALYSIS_PARMETIS,

  CS_PARAM_MUMPS_ANALYSIS_AUTO,

  CS_PARAM_MUMPS_N_ANALYSIS_ALGOS

} cs_param_mumps_analysis_algo_t;


/*! \enum cs_param_mumps_memory_usage_t
 *  \brief Strategy for the memory usage inside MUMPS
 *
 * \var CS_PARAM_MUMPS_MEMORY_CONSTRAINED
 * Strategy aiming at limiting the memory usage
 *
 * \var CS_PARAM_MUMPS_MEMORY_AUTO
 * Strategy relying on the default settings
 *
 * \var CS_PARAM_MUMPS_MEMORY_CPU_DRIVEN
 * Strategy aiming at the best CPU time
 */

typedef enum {

  CS_PARAM_MUMPS_MEMORY_CONSTRAINED,
  CS_PARAM_MUMPS_MEMORY_AUTO,
  CS_PARAM_MUMPS_MEMORY_CPU_DRIVEN,

  CS_PARAM_MUMPS_N_MEMORY_USAGES

} cs_param_mumps_memory_usage_t;


/*! \struct cs_param_mumps_t
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

  cs_param_mumps_analysis_algo_t   analysis_algo;
  cs_param_mumps_facto_type_t      facto_type;
  cs_param_mumps_memory_usage_t    mem_usage;

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

  bool    keep_ordering;    /*!< Mutualization of the ordering step */

  int     block_analysis;   /*!< Analysis is performed by block. Value of the
                                 block size. Not used if < 1 */

  int     ir_steps;         /*!< Number of steps for the Iterative Refinement */

} cs_param_mumps_t;

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create and initialize with the default settings a new structure
 *        storing a set of parameters used when calling MUMPS.
 *
 * \return a pointer to a new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_param_mumps_t *
cs_param_mumps_create(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy into a new structure the given set of parameters used when
 *        calling MUMPS
 *
 * \param[in] mumpsp   set of mumps parameters
 *
 * \return a pointer to a new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_param_mumps_t *
cs_param_mumps_copy(const cs_param_mumps_t  *mumpsp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log the structure storing the set of parameters used with MUMPS
 *
 * \param[in] name     name related to the current SLES
 * \param[in] mumpsp   set of mumps parameters
 */
/*----------------------------------------------------------------------------*/

void
cs_param_mumps_log(const char              *name,
                   const cs_param_mumps_t  *mumpsp);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PARAM_MUMPS_H__ */
