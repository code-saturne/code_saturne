/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2008 EDF S.A., France
 *
 *     contact: saturne-support@edf.fr
 *
 *     The Code_Saturne Kernel is free software; you can redistribute it
 *     and/or modify it under the terms of the GNU General Public License
 *     as published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version.
 *
 *     The Code_Saturne Kernel is distributed in the hope that it will be
 *     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with the Code_Saturne Kernel; if not, write to the
 *     Free Software Foundation, Inc.,
 *     51 Franklin St, Fifth Floor,
 *     Boston, MA  02110-1301  USA
 *
 *============================================================================*/

#ifndef __CS_SELECTOR_H__
#define __CS_SELECTOR_H__

/*============================================================================
 * Build selection lists for faces or cells
 *============================================================================*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Public functions definition for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Build a list of boundary faces verifying a given selection criteria.
 *----------------------------------------------------------------------------*/

void CS_PROCF(csgfbr, CSGFBR)
(
 const char   *const fstr,      /* <-- Fortran string */
 cs_int_t     *const len,       /* <-- String Length  */
 cs_int_t     *const n_faces,   /* --> number of faces */
 cs_int_t     *const face_list  /* --> face list  */
 CS_ARGF_SUPP_CHAINE
);

/*----------------------------------------------------------------------------
 * Build a list of interior faces verifying a given selection criteria.
 *----------------------------------------------------------------------------*/

void CS_PROCF(csgfac, CSGFAC)
(
 const char   *const fstr,      /* <-- Fortran string */
 cs_int_t     *const len,       /* <-- String Length  */
 cs_int_t     *const n_faces,   /* --> number of faces */
 cs_int_t     *const face_list  /* --> face list  */
 CS_ARGF_SUPP_CHAINE
);

/*----------------------------------------------------------------------------
 * Build a list of cells verifying a given selection criteria.
 *----------------------------------------------------------------------------*/

void CS_PROCF(csgcel, CSGCEL)
(
 const char   *const fstr,      /* <-- Fortran string */
 cs_int_t     *const len,       /* <-- String Length  */
 cs_int_t     *const n_cells,   /* --> number of cells */
 cs_int_t     *const cell_list  /* --> cell list  */
 CS_ARGF_SUPP_CHAINE
);

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Fill a list of boundary faces verifying a given selection criteria.
 *
 * parameters:
 *   criteria    <-- selection criteria string
 *   n_b_faces   --> number of selected interior faces
 *   b_face_list --> list of selected boundary faces
 *                   (1 to n, preallocated to cs_glob_mesh->n_b_faces)
 *----------------------------------------------------------------------------*/

void
cs_selector_get_b_face_list(const char  *criteria,
                            fvm_lnum_t  *n_b_faces,
                            fvm_lnum_t   b_face_list[]);

/*----------------------------------------------------------------------------
 * Fill a list of interior faces verifying a given selection criteria.
 *
 * parameters:
 *   criteria    <-- selection criteria string
 *   n_i_faces   --> number of selected interior faces
 *   i_face_list --> list of selected interior faces
 *                   (1 to n, preallocated to cs_glob_mesh->n_i_faces)
 *----------------------------------------------------------------------------*/

void
cs_selector_get_i_face_list(const char  *criteria,
                            fvm_lnum_t  *n_i_faces,
                            fvm_lnum_t   i_face_list[]);

/*----------------------------------------------------------------------------
 * Fill a list of cells verifying a given selection criteria.
 *
 * parameters:
 *   criteria  <-- selection criteria string
 *   n_cells   --> number of selected cells
 *   cell_list --> list of selected cells
 *                 (1 to n, preallocated to cs_glob_mesh->n_cells)
 *----------------------------------------------------------------------------*/

void
cs_selector_get_cell_list(const char  *criteria,
                          fvm_lnum_t  *n_cells,
                          fvm_lnum_t   cell_list[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SELECTOR_H__ */
