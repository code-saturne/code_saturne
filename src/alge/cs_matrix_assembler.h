#ifndef __CS_MATRIX_ASSEMBLER_H__
#define __CS_MATRIX_ASSEMBLER_H__

/*============================================================================
 * Incremental or general construction of matrix.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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

#include "cs_defs.h"

#include "cs_halo.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! Structure used to pre-build a matrix */

typedef struct _cs_matrix_assembler_t  cs_matrix_assembler_t;

/*============================================================================
 *  Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a matrix assembler structure.
 *
 * The associated matrix structure data will initially be empty, though
 * the range of rows associated with the current rank must be defined
 * immediately.
 *
 * \param[in]  l_range        global id range [min, max[ for local rank
 * \param[in]  separate_diag  if true, diagonal terms are handled separately
 *
 * \return  pointer to created matrix assembler structure;
 */
/*----------------------------------------------------------------------------*/

cs_matrix_assembler_t *
cs_matrix_assembler_create(cs_gnum_t  l_range[2],
                           bool       separate_diag);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy a matrix assembler structure.
 *
 * \param[in, out]  pointer to matrix assembler structure pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_assembler_destroy(cs_matrix_assembler_t  **ma);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add entries to a matrix assembler structure.
 *
 * This function should be called by a single thread for a given assembler.
 *
 * \param[in, out]  ma        pointer to matrix assembler structure
 * \param[in]       n         number of entries
 * \param[in]       g_col_id  global column ids associated with entries
 * \param[in]       g_row_id  global row ids associated with entries
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_assembler_add_ids(cs_matrix_assembler_t  *ma,
                            cs_lnum_t               n,
                            cs_gnum_t               g_row_id[],
                            cs_gnum_t               g_col_id[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute internal structures required by a matrix assembler.
 *
 * The associated vector halo is also computed at this stage.
 *
 * \param[in, out]  ma   pointer to matrix assembler structure
 *
 * This function should be called by a single thread for a given assembler.
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_assembler_compute_structure(cs_matrix_assembler_t  *ma);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign specific MPI communicator to matrix assembler.
 *
 * This must be called before \ref cs_matrix_assembler_compute.
 *
 * \param[in, out]  ma    pointer to matrix assembler structure
 * \param[in]       comm  assigned MPI communicator
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_MPI)

void
cs_matrix_assembler_set_comm(cs_matrix_assembler_t  *ma,
                             MPI_Comm                comm);

#endif /* HAVE_MPI */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a pointer to the halo structure associated with a
 *        matrix assembler.
 *
 * \param[in]  ma   pointer to matrix assembler structure
 *
 * \return  pointer to halo structure
 */
/*----------------------------------------------------------------------------*/

const cs_halo_t *
cs_matrix_assembler_get_halo(const cs_matrix_assembler_t  *ma);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MATRIX_ASSEMBLER_H__ */
