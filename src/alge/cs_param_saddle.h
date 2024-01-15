#ifndef __CS_PARAM_SADDLE_H__
#define __CS_PARAM_SADDLE_H__

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

#include "cs_param_sles.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
 * \file cs_param_saddle.h

 * \brief Handle the settings of saddle-point systems.
 *        These systems arise from the monolithic coupling of the Navier-Stokes
 *        equations or in mixed formulation of scalar-valued equations.
 */

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*!
 * \enum cs_param_saddle_precond_t
 *
 * \brief Type of preconditioner used to solve a saddle-point system. Up to
 * now, this happens only with CDO cell-based schemes. Saddle-point system
 * arising from the Stokes or Navier-Stokes equations are handled differently
 * even though there are several similarities.
 *
 * \var CS_PARAM_SADDLE_PRECOND_NONE
 * No preconditioner to apply.
 *
 * \var CS_PARAM_SADDLE_PRECOND_DIAG_SCHUR
 * A block-diagonal preconditioner is used. The (1,1)-block is an approximation
 * of the (1,1)-block in the saddle-point system based on a cheap
 * resolution. The parameter settings for this resolution relies on the
 * structure \ref cs_param_sles_t. The (2,2)-block in the preconditioner is the
 * Schur approximation. Please refer to \ref cs_param_schur_approx_t to get a
 * detailed view of the possibilities.
 *
 * \var CS_PARAM_SADDLE_PRECOND_LOWER_SCHUR
 * A 2x2 block matrix is used as preconditioner with the (1,2)-block fills with
 * zero. The (1,1)-block is an approximation of the (1,1)-block in the
 * saddle-point system based on a cheap resolution. The parameter settings for
 * this resolution relies on the structure \ref cs_param_sles_t The (2,2)-block
 * in the preconditioner is the Schur approximation. Please refer to \ref
 * cs_param_schur_approx_t to get a detailed view of the possibilities.
 *
 * \var CS_PARAM_SADDLE_PRECOND_UPPER_SCHUR
 * A 2x2 block matrix is used as preconditioner with the (2,1)-block fills with
 * zero. The (1,1)-block is an approximation of the (1,1)-block in the
 * saddle-point system based on a cheap resolution. The parameter settings for
 * this resolution relies on the structure \ref cs_param_sles_t The (2,2)-block
 * in the preconditioner is the Schur approximation. Please refer to \ref
 * cs_param_schur_approx_t to get a detailed view of the possibilities.
 */

typedef enum {

  CS_PARAM_SADDLE_PRECOND_NONE,
  CS_PARAM_SADDLE_PRECOND_DIAG_SCHUR,
  CS_PARAM_SADDLE_PRECOND_LOWER_SCHUR,
  CS_PARAM_SADDLE_PRECOND_UPPER_SCHUR,

  CS_PARAM_SADDLE_N_PRECOND

} cs_param_saddle_precond_t;

/*!
 * \enum cs_param_saddle_solver_t
 *
 * \brief Type of solver used to solve a saddle-point system. Up to now, this
 * happens only with CDO cell-based schemes. Saddle-point system arising from
 * the Stokes or Navier-Stokes equations are handled differently even though
 * there are several similarities. The main differences are the fact that the
 * (1,1) block is vector-valued and that there can be specific optimizations in
 * the preconditioning (the Schur complement approximation for instance) since
 * one has a more advanced knowledge of the system.
 *
 * \var CS_PARAM_SADDLE_SOLVER_NONE
 * No solver defined. No saddle-point to solve.
 *
 * \var CS_PARAM_SADDLE_SOLVER_GCR
 * Iterative solver for indefinite systems. This solver relies on a specific
 * storage of the saddle-point system: (1,1)-block is assembled and the
 * (2,1)-block is unassembled. The (1,2) is not stored since this the
 * transposition of the (2,1)-block.
 *
 * \var CS_PARAM_SADDLE_SOLVER_MINRES
 * Iterative solver for indefinite symmetric systems. (This is like a CG solver
 * for symmetric saddle-point systems). This solver relies on a specific
 * storage of the saddle-point system: (1,1)-block is assembled and the
 * (2,1)-block is unassembled. The (1,2) is not stored since this the
 * transposition of the (2,1)-block.
 *
 * \var CS_PARAM_SADDLE_SOLVER_MUMPS
 * Sparse direct solver based on the external library called MUMPS
 */

typedef enum {

  CS_PARAM_SADDLE_SOLVER_NONE,

  CS_PARAM_SADDLE_SOLVER_GCR,
  CS_PARAM_SADDLE_SOLVER_MINRES,
  CS_PARAM_SADDLE_SOLVER_MUMPS,

  CS_PARAM_SADDLE_N_SOLVERS

} cs_param_saddle_solver_t;

/*!
 * \struct cs_param_saddle_t
 * \brief Structure storing all metadata related to the resolution of a
 *        saddle-point linear system. A saddle-point system is depicted as
 *
 * |  A  :  Bt |  where A is the (1,1) block
 * |-----:-----|        Bt is the (1,2) block which the transposed operator
 * |  B  :  0  |        w.r.t. B, the (2,1)-block
 *
 */

typedef struct {

  /*!
    * \var verbosity
    * verbosity (level of information displayed)
    */

  int                         verbosity;

  /*!
    * \var name
    * name of the saddle system or NULL (not mandatory). If NULL, the name
    * of the (1,1)-block is used.
    */

  char                       *name;

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

  /*! \var block11_sles_param
   * Set of parameters used to solve (1,1)-block i.e. the A matrix. This
   * is shared with the \ref cs_equation_param_t structure.
   */

  const cs_param_sles_t      *block11_sles_param;

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
   * Set of parameters used to solve the Schur complement if needed. This
   * depends on the type of Schur approximation which has been chosen.
   */

  cs_param_sles_t            *schur_sles_param;

  /*! @} */

} cs_param_saddle_t;

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a cs_param_saddle_t structure and assign default settings
 *
 * \param[in] block11_slesp   set of parameters for the (1,1) block
 *
 * \return a pointer to the new cs_param_saddle_t structure
 */
/*----------------------------------------------------------------------------*/

cs_param_saddle_t *
cs_param_saddle_create(const cs_param_sles_t   *block11_slesp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free the structure storing the parameter settings for a saddle-point
 *        system
 *
 * \param[in, out] p_saddlep    double pointer to the structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_param_saddle_free(cs_param_saddle_t  **p_saddlep);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the name of the saddle-point system.
 *
 * \param[in]      basename   prefix for the naming of the Schur system
 * \param[in, out] saddlep    pointer to the structure to update
 */
/*----------------------------------------------------------------------------*/

void
cs_param_saddle_set_name(const char         *name,
                         cs_param_saddle_t  *saddlep);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize a \ref cs_param_sles_t structure for the Schur
 *        approximation nested inside a \ref cs_param_saddle_t structure. By
 *        default, this member is not allocated. Do nothing if the related
 *        structure is already allocated.
 *
 * \param[in, out] saddlep    pointer to the structure to update
 */
/*----------------------------------------------------------------------------*/

void
cs_param_saddle_init_schur_sles(cs_param_saddle_t  *saddlep);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy a cs_param_saddle_t structure from ref to dest
 *
 * \param[in]      ref     reference structure to be copied
 * \param[in, out] dest    destination structure
 */
/*----------------------------------------------------------------------------*/

void
cs_param_saddle_copy(const cs_param_saddle_t  *ref,
                     cs_param_saddle_t        *dest);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log the setup information for the given cs_param_saddle_t structure
 *
 * \param[in] saddlep     pointer to the structure
 */
/*----------------------------------------------------------------------------*/

void
cs_param_saddle_log(const cs_param_saddle_t  *saddlep);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PARAM_SADDLE_H__ */
