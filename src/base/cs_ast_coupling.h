#ifndef __CS_AST_COUPLING_H__
#define __CS_AST_COUPLING_H__

/*============================================================================
 * code_aster coupling
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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*============================================================================
 * Structure definition
 *============================================================================*/

typedef struct _cs_ast_coupling_t  cs_ast_coupling_t;

/*============================================================================
 * Global variable definitions
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initial exchange with code_aster
 *
 * \param[in]  verbosity      verbosity level for code_aster coupling
 * \param[in]  visualization  visualization level for code_aster coupling
 * \param[in]  nalimx         maximum number of implicitation iterations of
 *                            the structure displacement
 * \param[in]  epalim         relative precision of implicitation of
 *                            the structure displacement
 */
/*----------------------------------------------------------------------------*/

void
cs_ast_coupling_initialize(int        verbosity,
                           int        visualization,
                           int        nalimx,
                           cs_real_t  epalim);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize exchange with code_aster
 */
/*----------------------------------------------------------------------------*/

void
cs_ast_coupling_finalize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Extract and exchange mesh information for surfaces coupled with
 *        code_aster.
 *
 * \param[in]  n_faces   number of coupled faces.
 * \param[in]  face_ids  ids of coupled faces (ordered by increasing id)
 * \param[in]  almax     characteristic macroscopic domain length
 */
/*----------------------------------------------------------------------------*/

void
cs_ast_coupling_geometry(cs_lnum_t         n_faces,
                         const cs_lnum_t  *face_ids,
                         cs_real_t         almax);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Exchange time-step information with code_aster.
 *
 * \param[in, out]  c_dt  time step at each cell
 */
/*----------------------------------------------------------------------------*/

void
cs_ast_coupling_exchange_time_step(cs_real_t  c_dt[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Send stresses acting on the fluid/structure interface
 *        and receive displacements.
 *
 * \param[in]  fluid_forces  forces from fluid at coupled faces
 */
/*----------------------------------------------------------------------------*/

void
cs_ast_coupling_exchange_fields(const cs_real_t  fluid_forces[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute predicted or exact displacement of the
 *        fluid/structure interface.
 *
 * \param[out]  disp  prescribed displacement at vertices
 */
/*----------------------------------------------------------------------------*/

void
cs_ast_coupling_compute_displacement(cs_real_t  disp[][3]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Receive convergence value of code_saturne/code_aster coupling
 *
 * \return  convergence indicator computed by coupling scheme
 *          (1: converged, 0: not converged)
 */
/*----------------------------------------------------------------------------*/

int
cs_ast_coupling_get_ext_cvg(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Send global convergence value of FSI calculations
 *
 * \param[in]  icved  convergence indicator (1: converged, 0: not converged)
 */
/*----------------------------------------------------------------------------*/

void
cs_ast_coupling_send_cvg(int  icved);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_AST_COUPLING_H__ */
