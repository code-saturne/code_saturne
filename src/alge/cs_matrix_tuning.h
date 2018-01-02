#ifndef __CS_MATRIX_TUNING_H__
#define __CS_MATRIX_TUNING_H__

/*============================================================================
 * Sparse Matrix Representation and Operations
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Tune local matrix.vector product operations.
 *
 * To avoid multiplying structures for multiple matrix fill-ins,
 * an array of tuning types may be provided, and weights may be
 * associated to each type based on the expected usage of each fill-in
 * type. If n_fill_types is set to 0, these arrays are ignored, and their
 * following default is used:
 *
 *   CS_MATRIX_SCALAR      0.5
 *   CS_MATRIX_SCALAR_SYM  0.25
 *   CS_MATRIX_33_BLOCK_D  0.25
 *
 * parameters:
 *   t_measure      <-- minimum time for each measure
 *   n_types        <-- number of matrix types tuned for, or 0
 *   n_fill_types   <-- number of fill types tuned for, or 0
 *   types          <-- array of matrix types tuned for, or NULL
 *   fill_types     <-- array of fill types tuned for, or NULL
 *   fill_weights   <-- weight of fill types tuned for, or NULL
 *   n_min_products <-- minimum number of SpMv products (to estimate
 *                      amortization of coefficients assignment)
 *   n_cells        <-- number of local cells
 *   n_cells_ext    <-- number of cells including ghost cells (array size)
 *   n_faces        <-- local number of internal faces
 *   face_cell      <-- face -> cells connectivity
 *   halo           <-- cell halo structure
 *   numbering      <-- vectorization or thread-related numbering info, or NULL
 *
 * returns:
 *   pointer to tuning results structure
 *----------------------------------------------------------------------------*/

cs_matrix_variant_t *
cs_matrix_variant_tuned(double                 t_measure,
                        int                    n_types,
                        int                    n_fill_types,
                        cs_matrix_type_t       types[],
                        cs_matrix_fill_type_t  fill_types[],
                        double                 fill_weights[],
                        int                    n_min_products,
                        cs_lnum_t              n_cells,
                        cs_lnum_t              n_cells_ext,
                        cs_lnum_t              n_faces,
                        const cs_lnum_2_t     *face_cell,
                        const cs_halo_t       *halo,
                        const cs_numbering_t  *numbering);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MATRIX_TUNING_H__ */
