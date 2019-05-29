#ifndef __CS_EQUATION_ASSEMBLE_H__
#define __CS_EQUATION_ASSEMBLE_H__

/*============================================================================
 * Routines to handle the assembly process of equatino discretized with CDO
 * schemes
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_cdo_connect.h"
#include "cs_cdo_local.h"
#include "cs_matrix.h"
#include "cs_matrix_assembler.h"
#include "cs_param.h"
#include "cs_param_cdo.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _cs_equation_assemble_t  cs_equation_assemble_t;

/*============================================================================
 * Function pointer type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise system into the global algebraic system
 *         Block or no block versions are handled
 *
 * \param[in]      csys     cellwise view of the algebraic system
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] eqa      pointer to an equation assembly structure
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_equation_assembly_t)(const cs_cell_sys_t                    *csys,
                         const cs_range_set_t                   *rset,
                         cs_equation_assemble_t                 *eqa,
                         cs_matrix_assembler_values_t           *mav);

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the pointer to a requested \ref cs_matrix_structure_t
 *         structure
 *
 * \param[in]  flag_id       id in the array of matrix structures
 *
 * \return  a pointer to a cs_matrix_structure_t
 */
/*----------------------------------------------------------------------------*/

cs_matrix_structure_t *
cs_equation_get_matrix_structure(int  flag);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get a pointer to a cs_equation_assemble_t structure related
 *         to a given thread
 *
 * \param[in]  t_id    id in the array of pointer
 *
 * \return a pointer to a cs_equation_assemble_t structure
 */
/*----------------------------------------------------------------------------*/

cs_equation_assemble_t *
cs_equation_assemble_get(int    t_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize matrix-related structures according to
 *         the type of discretization used for this simulation
 *
 * \param[in]  connect      pointer to a cs_cdo_connect_t structure
 * \param[in]  vb_flag      metadata for Vb schemes
 * \param[in]  vcb_flag     metadata for V+C schemes
 * \param[in]  fb_flag      metadata for Fb schemes
 * \param[in]  hho_flag     metadata for HHO schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_assemble_init(const cs_cdo_connect_t       *connect,
                          cs_flag_t                     vb_flag,
                          cs_flag_t                     vcb_flag,
                          cs_flag_t                     fb_flag,
                          cs_flag_t                     hho_flag);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free matrix-related structures used during the simulation.
 *         Display overall statistic about the assembly stage for CDO schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_assemble_finalize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the function pointer used to assemble the algebraic system
 *
 * \param[in] scheme     space discretization scheme
 * \param[in] ma_id      id in the array of matrix assembler
 *
 * \return a function pointer cs_equation_assembly_t
 */
/*----------------------------------------------------------------------------*/

cs_equation_assembly_t *
cs_equation_assemble_set(cs_param_space_scheme_t    scheme,
                         int                        ma_id);

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise system into the global algebraic system
 *         Scalar-valued case. Parallel and with openMP threading.
 *
 * \param[in]      csys     cellwise view of the algebraic system
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] eqa      pointer to a matrix assembler buffers
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_assemble_matrix_mpit(const cs_cell_sys_t              *csys,
                                 const cs_range_set_t             *rset,
                                 cs_equation_assemble_t           *eqa,
                                 cs_matrix_assembler_values_t     *mav);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise system into the global algebraic system
 *         Scalar-valued case. Parallel without openMP threading.
 *
 * \param[in]      csys     cellwise view of the algebraic system
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] eqa      pointer to a matrix assembler buffers
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_assemble_matrix_mpis(const cs_cell_sys_t              *csys,
                                 const cs_range_set_t             *rset,
                                 cs_equation_assemble_t           *eqa,
                                 cs_matrix_assembler_values_t     *mav);

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise system into the global algebraic system
 *         Scalar-valued case. Sequential and with openMP threading.
 *
 * \param[in]      csys     cellwise view of the algebraic system
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] eqa      pointer to a matrix assembler buffers
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_assemble_matrix_seqt(const cs_cell_sys_t             *csys,
                                 const cs_range_set_t            *rset,
                                 cs_equation_assemble_t          *eqa,
                                 cs_matrix_assembler_values_t    *mav);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise system into the global algebraic system
 *         Scalar-valued case. Sequential and without openMP threading.
 *
 * \param[in]      csys     cellwise view of the algebraic system
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] eqa      pointer to a matrix assembler buffers
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_assemble_matrix_seqs(const cs_cell_sys_t              *csys,
                                 const cs_range_set_t             *rset,
                                 cs_equation_assemble_t           *eqa,
                                 cs_matrix_assembler_values_t     *mav);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise system into the global algebraic system.
 *         Case of a block 3x3 entries. Expand each row.
 *         Sequential run without openMP threading.
 *
 * \param[in]      csys     cellwise view of the algebraic system
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] eqa      pointer to an equation assembly structure
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_assemble_eblock33_matrix_seqs(const cs_cell_sys_t           *csys,
                                          const cs_range_set_t          *rset,
                                          cs_equation_assemble_t        *eqa,
                                          cs_matrix_assembler_values_t  *mav);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise system into the global algebraic system.
 *         Case of a block 3x3 entries. Expand each row.
 *         Sequential run with openMP threading.
 *
 * \param[in]      csys     cellwise view of the algebraic system
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] eqa      pointer to an equation assembly structure
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_assemble_eblock33_matrix_seqt(const cs_cell_sys_t           *csys,
                                          const cs_range_set_t          *rset,
                                          cs_equation_assemble_t        *eqa,
                                          cs_matrix_assembler_values_t  *mav);

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise system into the global algebraic system.
 *         Case of a block 3x3 entries. Expand each row.
 *         Parallel run without openMP threading.
 *
 * \param[in]      csys     cellwise view of the algebraic system
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] eqa      pointer to an equation assembly structure
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_assemble_eblock33_matrix_mpis(const cs_cell_sys_t           *csys,
                                          const cs_range_set_t          *rset,
                                          cs_equation_assemble_t        *eqa,
                                          cs_matrix_assembler_values_t  *mav);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise system into the global algebraic system.
 *         Case of a block 3x3 entries. Expand each row.
 *         Parallel run with openMP threading.
 *
 * \param[in]      csys     cellwise view of the algebraic system
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] eqa      pointer to an equation assembly structure
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_assemble_eblock33_matrix_mpit(const cs_cell_sys_t           *csys,
                                          const cs_range_set_t          *rset,
                                          cs_equation_assemble_t        *eqa,
                                          cs_matrix_assembler_values_t  *mav);

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise system into the global algebraic system.
 *         Case of a block NxN entries. Expand each row.
 *         Sequential run without openMP threading.
 *
 * \param[in]      csys     cellwise view of the algebraic system
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] eqa      pointer to an equation assembly structure
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_assemble_eblock_matrix_seqs(const cs_cell_sys_t           *csys,
                                        const cs_range_set_t          *rset,
                                        cs_equation_assemble_t        *eqa,
                                        cs_matrix_assembler_values_t  *mav);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise system into the global algebraic system.
 *         Case of a block NxN entries. Expand each row.
 *         Sequential run with openMP threading.
 *
 * \param[in]      csys     cellwise view of the algebraic system
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] eqa      pointer to an equation assembly structure
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_assemble_eblock_matrix_seqt(const cs_cell_sys_t           *csys,
                                        const cs_range_set_t          *rset,
                                        cs_equation_assemble_t        *eqa,
                                        cs_matrix_assembler_values_t  *mav);

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise system into the global algebraic system.
 *         Case of a block NxN entries. Expand each row.
 *         Parallel run without openMP threading.
 *
 * \param[in]      csys     cellwise view of the algebraic system
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] eqa      pointer to an equation assembly structure
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_assemble_eblock_matrix_mpis(const cs_cell_sys_t           *csys,
                                        const cs_range_set_t          *rset,
                                        cs_equation_assemble_t        *eqa,
                                        cs_matrix_assembler_values_t  *mav);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise system into the global algebraic system.
 *         Case of a block NxN entries. Expand each row.
 *         Parallel run with openMP threading.
 *
 * \param[in]      csys     cellwise view of the algebraic system
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] eqa      pointer to an equation assembly structure
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_assemble_eblock_matrix_mpit(const cs_cell_sys_t           *csys,
                                        const cs_range_set_t          *rset,
                                        cs_equation_assemble_t        *eqa,
                                        cs_matrix_assembler_values_t  *mav);

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_EQUATION_ASSEMBLE_H__ */
