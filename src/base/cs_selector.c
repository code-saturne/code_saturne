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

#if defined(HAVE_CONFIG_H)
#include "cs_config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <bft_mem_usage.h>
#include <bft_mem.h>
#include <bft_error.h>
#include <bft_printf.h>
#include <bft_sys_info.h>
#include <bft_timer.h>

#include "cs_mesh.h"
#include "cs_selector.h"

#include "fvm_selector.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

#define CS_SELECTOR_STR_LEN 127

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 *  Public function definitions for Fortran API
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
)
{
  char _c_string[CS_SELECTOR_STR_LEN + 1];
  char *c_string = _c_string;
  int i;
  int c_len = *len -1;

  /* Initialization */

  *n_faces = 0;

  /* Copy fstr without last blanks  */
  while(fstr[c_len--] == ' ' &&  c_len >= 0);

  if (c_len < -1)
    return;

  c_len += 2;

  if (c_len > CS_SELECTOR_STR_LEN)
    BFT_MALLOC(c_string, c_len + 1, char);

  for(i = 0; i < c_len; i++)
    c_string[i] = fstr[i];
  c_string[c_len] = '\0';

  /* Get faces with C string */

  cs_selector_get_b_face_list(c_string, n_faces, face_list);

  if (c_string != _c_string)
    BFT_FREE(c_string);
}

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
)
{
  char _c_string[CS_SELECTOR_STR_LEN + 1];
  char *c_string = _c_string;
  int i;
  int c_len = *len -1;

  /* Initialization */

  *n_faces = 0;

  /* Copy fstr without last blanks  */
  while(fstr[c_len--] == ' ' &&  c_len >= 0);

  if (c_len < -1)
    return;

  c_len += 2;

  if (c_len > CS_SELECTOR_STR_LEN)
    BFT_MALLOC(c_string, c_len + 1, char);

  for(i = 0; i < c_len; i++)
    c_string[i] = fstr[i];
  c_string[c_len] = '\0';

  /* Get faces with C string */

  cs_selector_get_i_face_list(c_string, n_faces, face_list);

  if (c_string != _c_string)
    BFT_FREE(c_string);
}

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
)
{
  char _c_string[CS_SELECTOR_STR_LEN + 1];
  char *c_string = _c_string;
  int i;
  int c_len = *len -1;

  /* Initialization */

  *n_cells = 0;

  /* Copy fstr without last blanks  */
  while(fstr[c_len--] == ' ' &&  c_len >= 0);

  if (c_len < -1)
    return;

  c_len += 2;

  if (c_len > CS_SELECTOR_STR_LEN)
    BFT_MALLOC(c_string, c_len + 1, char);

  for(i = 0; i < c_len; i++)
    c_string[i] = fstr[i];
  c_string[c_len] = '\0';

  /* Get cells with C string */

  cs_selector_get_cell_list(c_string, n_cells, cell_list);

  if (c_string != _c_string)
    BFT_FREE(c_string);
}

/*=============================================================================
 * Public function definitions
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
                            fvm_lnum_t   b_face_list[])
{
  int c_id;

  *n_b_faces = 0;

  c_id = fvm_selector_get_list(cs_glob_mesh->select_b_faces,
                               criteria,
                               n_b_faces,
                               b_face_list);

  if (fvm_selector_n_missing(cs_glob_mesh->select_b_faces, c_id) > 0) {
    const char *missing
      = fvm_selector_get_missing(cs_glob_mesh->select_b_faces, c_id, 0);
    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_("The group or attribute \"%s\" in the selection\n"
                 "criteria:\n"
                 "\"%s\"\n"
                 " does not correspond to any boundary face.\n"),
               missing, criteria);
  }
}

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
                            fvm_lnum_t   i_face_list[])
{
  int c_id;

  *n_i_faces = 0;

  c_id = fvm_selector_get_list(cs_glob_mesh->select_i_faces,
                               criteria,
                               n_i_faces,
                               i_face_list);

  if (fvm_selector_n_missing(cs_glob_mesh->select_i_faces, c_id) > 0) {
    const char *missing
      = fvm_selector_get_missing(cs_glob_mesh->select_i_faces, c_id, 0);
    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_("The group or attribute \"%s\" in the selection\n"
                 "criteria:\n"
                 "\"%s\"\n"
                 " does not correspond to any interior face.\n"),
               missing, criteria);
  }
}

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
                          fvm_lnum_t   cell_list[])
{
  int c_id;

  *n_cells = 0;

  c_id = fvm_selector_get_list(cs_glob_mesh->select_cells,
                               criteria,
                               n_cells,
                               cell_list);

  if (fvm_selector_n_missing(cs_glob_mesh->select_cells, c_id) > 0) {
    const char *missing
      = fvm_selector_get_missing(cs_glob_mesh->select_cells, c_id, 0);
    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_("The group or attribute \"%s\" in the selection\n"
                 "criteria:\n"
                 "\"%s\"\n"
                 " does not correspond to any cell.\n"),
               missing, criteria);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
