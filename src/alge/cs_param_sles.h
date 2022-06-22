#ifndef __CS_PARAM_SLES_H__
#define __CS_PARAM_SLES_H__

/*============================================================================
 * Routines to handle the SLES settings
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

/*!
 * \struct cs_param_sles_t
 * \brief Structure storing all metadata related to the resolution of a linear
 *        system with an iterative solver.
 */

typedef struct {

  bool                     setup_done;   /*!< SLES setup step has been done */
  int                      verbosity;    /*!< SLES verbosity */
  int                      field_id;     /*!< Field id related to a SLES
                                           By default, this is set to -1 */
  char                    *name;         /*!< System name (equation name if
                                           this is automatic) */

  cs_param_sles_class_t      solver_class; /*!< class of SLES to consider  */
  cs_param_precond_type_t    precond;      /*!< type of preconditioner */
  cs_param_itsol_type_t      solver;       /*!< type of solver */
  bool                       flexible;     /*!< need a flexible variant ? */

  /*! \var amg_type
   *  type of AMG algorithm if needed */

  cs_param_amg_type_t        amg_type;

  /*! \var pcd_block_type
   *  type of block preconditioner to use (only meaningful for vector-valued
   *  systems or more complex systems */

  cs_param_precond_block_t   pcd_block_type;

  /*! \var resnorm_type
   *  normalized or not the norm of the residual used for the stopping criterion
   *  See \ref CS_EQKEY_ITSOL_RESNORM_TYPE for more details. */

  cs_param_resnorm_type_t  resnorm_type;
  int                      n_max_iter;   /*!< max. number of iterations */
  int                      restart;      /*!< max. iter. before restarting  */
  double                   eps;          /*!< stopping criterion on accuracy */

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
