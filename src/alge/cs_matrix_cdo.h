#ifndef __CS_MATRIX_CDO_H__
#define __CS_MATRIX_CDO_H__

/*============================================================================
 * Sparse Matrix Representation and Operations
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2015 EDF S.A.

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
#include "cs_numbering.h"
#include "cs_halo_perio.h"
#include "cs_matrix_priv.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Structure associated with opaque matrix structure object */

typedef struct _cs_matrix_cdo_structure_t cs_matrix_cdo_structure_t;

/*============================================================================
 *  Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * \brief Create a CDO matrix structure.
 *
 * Only CSR and MSR formats are handled.
 * col_id is sorted row by row during the creation of this structure.
 *
 * \param[in] type       type of matrix considered
 * \param[in] owner      deallocate members in structure only if owner=true
 * \param[in] have_diag  indicates if the diagonal structure contains nonzeroes
 * \param[in] n_rows     local number of rows
 * \param[in] n_rows_ext local number of rows + ghost rows implied in sync.
 * \param[in] n_cols     local number of cols
 * \param[in] idx        index on rows
 * \param[inout] col_id  array of colum ids related to the row index
 * \param[in] halo       halo structure for synchronization in parallel, or NULL
 * \param[in] numbering  vectorization or thread-related numbering info, or NULL
 *
 * \return  a pointer to a created CDO matrix structure
 *----------------------------------------------------------------------------*/

cs_matrix_cdo_structure_t *
cs_matrix_cdo_structure_create(cs_matrix_type_t       type,
                               bool                   owner,
                               bool                   have_diag,
                               cs_lnum_t              n_rows,
                               cs_lnum_t              n_rows_ext,
                               cs_lnum_t              n_cols,
                               cs_lnum_t             *idx,
                               cs_lnum_t             *col_id,
                               const cs_halo_t       *halo,
                               const cs_numbering_t  *numbering);

/*----------------------------------------------------------------------------
 * \brief Destroy a CDO matrix structure.
 *
 * \param[inout]  ms  pointer to a CDO matrix structure pointer
 *----------------------------------------------------------------------------*/

void
cs_matrix_cdo_structure_destroy(cs_matrix_cdo_structure_t  **ms);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a matrix container using a given CDO matrix structure.
 *
 * Note that the matrix container maps to the assigned structure,
 * so it must be destroyed before that structure.
 *
 * \param[in]  ms  associated CDO matrix structure
 *
 * \return  pointer to created matrix structure;
 */
/*----------------------------------------------------------------------------*/

cs_matrix_t *
cs_matrix_cdo_create(const cs_matrix_cdo_structure_t  *ms);


/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MATRIX_CDO_H__ */
