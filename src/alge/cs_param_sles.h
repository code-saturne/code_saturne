#ifndef __CS_PARAM_SLES_H__
#define __CS_PARAM_SLES_H__

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
  switch (solver) {

  case CS_PARAM_ITSOL_MUMPS:
  case CS_PARAM_ITSOL_MUMPS_FLOAT:
  case CS_PARAM_ITSOL_MUMPS_FLOAT_LDLT:
  case CS_PARAM_ITSOL_MUMPS_FLOAT_SYM:
  case CS_PARAM_ITSOL_MUMPS_LDLT:
  case CS_PARAM_ITSOL_MUMPS_SYM:
    return true;

  default:
    return false;

  }
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
 * \brief  Free a \ref cs_param_sles_t structure
 *
 * \param[in, out]  slesp    pointer to a \cs_param_sles_t structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_free(cs_param_sles_t   **p_slesp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Log information related to the linear settings stored in the
 *         structure
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
 * \param[in]       src    reference cs_param_sles_t structure to copy
 * \param[in, out]  dst    copy of the reference at exit
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_copy_from(cs_param_sles_t   *src,
                        cs_param_sles_t   *dst);

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
 * \brief  Check if the setting related to the AMG is consistent with the
 *         solver class.
 *
 * \param[in, out] slesp    pointer to a cs_pparam_sles_t structure
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
