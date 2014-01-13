#ifndef __CS_SELECTOR_H__
#define __CS_SELECTOR_H__

/*============================================================================
 * Build selection lists for faces or cells
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2014 EDF S.A.

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

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Public function prototypes for Fortran API
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

/*----------------------------------------------------------------------------
 * Build a list of cells verifying a given selection criteria.
 *----------------------------------------------------------------------------*/

void CS_PROCF(csgceb, CSGCEB)
(
 const char   *const fstr,        /* <-- Fortran string */
 cs_int_t     *const len,         /* <-- String Length  */
 cs_int_t     *const n_i_faces,   /* --> number of interior faces */
 cs_int_t     *const n_b_faces,   /* --> number of boundary faces */
 cs_int_t     *const i_face_list, /* --> interior face list  */
 cs_int_t     *const b_face_list  /* --> boundary face list  */
 CS_ARGF_SUPP_CHAINE
);

/*----------------------------------------------------------------------------
 * Build a list of interior faces belonging to a given periodicity.
 *----------------------------------------------------------------------------*/

void CS_PROCF(getfpe, GETFPE)
(
 cs_int_t     *const perio_num, /* <-- Periodicity number */
 cs_int_t     *const n_faces,   /* --> number of faces */
 cs_int_t     *const face_list  /* --> face list  */
 CS_ARGF_SUPP_CHAINE
);

/*----------------------------------------------------------------------------
 * Build a list of families verifying a given selection criteria.
 *----------------------------------------------------------------------------*/

void CS_PROCF(csgfam, CSGFAM)
(
 const char   *const fstr,         /* <-- Fortran string */
 cs_int_t     *const len,          /* <-- String Length  */
 cs_int_t     *const n_families,   /* --> number of families */
 cs_int_t     *const family_list   /* --> family list  */
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
                            cs_lnum_t   *n_b_faces,
                            cs_lnum_t    b_face_list[]);

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
                            cs_lnum_t   *n_i_faces,
                            cs_lnum_t    i_face_list[]);

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
                          cs_lnum_t   *n_cells,
                          cs_lnum_t    cell_list[]);

/*----------------------------------------------------------------------------
 * Fill lists of faces at the boundary of a set of cells verifying a given
 * selection criteria.
 *
 * parameters:
 *   criteria    <-- selection criteria string
 *   n_i_faces   --> number of selected interior faces
 *   n_b_faces   --> number of selected interior faces
 *   i_face_list --> list of selected interior faces
 *                   (1 to n, preallocated to cs_glob_mesh->n_i_faces)
 *   b_face_list --> list of selected boundary faces
 *                   (1 to n, preallocated to cs_glob_mesh->n_b_faces)
 *----------------------------------------------------------------------------*/

void
cs_selector_get_cells_boundary(const char  *criteria,
                               cs_lnum_t   *n_i_faces,
                               cs_lnum_t   *n_b_faces,
                               cs_lnum_t    i_face_list[],
                               cs_lnum_t    b_face_list[]);

/*----------------------------------------------------------------------------
 * Fill a list of interior faces belonging to a given periodicity.
 *
 * parameters:
 *   perio_num   <-- periodicity number
 *   n_i_faces   --> number of selected interior faces
 *   i_face_list --> list of selected interior faces
 *                   (1 to n, preallocated to cs_glob_mesh->n_i_faces)
 *----------------------------------------------------------------------------*/

void
cs_selector_get_perio_face_list(int          perio_num,
                                cs_lnum_t   *n_i_faces,
                                cs_lnum_t    i_face_list[]);

/*----------------------------------------------------------------------------
 * Fill a list of families verifying a given selection criteria.
 *
 * parameters:
 *   criteria    <-- selection criteria string
 *   n_families  --> number of selected families
 *   family_list --> list of selected families
 *                   (0 to n, preallocated to cs_glob_mesh->n_families + 1)
 *----------------------------------------------------------------------------*/

void
cs_selector_get_family_list(const char  *criteria,
                            cs_lnum_t   *n_families,
                            cs_int_t     family_list[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SELECTOR_H__ */
