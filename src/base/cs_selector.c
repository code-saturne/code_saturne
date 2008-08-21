/*============================================================================
 *
 *                    Code_Saturne version 1.3
 *                    ------------------------
 *
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

/* pour les API Fortran : a deplacer */
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_selector.h"

#include "fvm_selector.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/* API Fortran : a deplacer */

/*----------------------------------------------------------------------------
 * Give the list of the faces checking the constraint defined by
 * the Fortran string "fstr"
 *----------------------------------------------------------------------------*/

void CS_PROCF(csgfbr, CSGFBR)
(
 const char          *const fstr,     /* <-- Fortran string */
 int                 *const len,      /* <-- String Length  */
 int                 *const n_faces,  /* --> number of faces */
 int                 *const faces     /* --> faces */
 CS_ARGF_SUPP_CHAINE
)
{
  char * cpyfstr;
  int i, c_id;
  int lenwithoutblank;

  /* Initialization */
  *n_faces = 0;

  /* Copy fstr without last blanks  */
  lenwithoutblank = *len - 1;
  while(fstr[lenwithoutblank--] == ' ' &&  lenwithoutblank >= 0);

  if (lenwithoutblank < -1)
    return;

  lenwithoutblank += 2;

  BFT_MALLOC(cpyfstr, lenwithoutblank + 1, char);

  for(i = 0 ; i < lenwithoutblank; i++)
    cpyfstr[i] = fstr[i];
  cpyfstr[lenwithoutblank] = '\0';

  /* Get faces with C string */

  c_id = fvm_selector_get_list(cs_glob_mesh->select_b_faces,
                               cpyfstr,
                               n_faces,
                               faces);

  if (fvm_selector_n_missing(cs_glob_mesh->select_b_faces, c_id) > 0) {
    const char *missing
      = fvm_selector_get_missing(cs_glob_mesh->select_b_faces, c_id, 0);
    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_("The group or attribute \"%s\" in the selection\n"
                   "criteria:\n"
                   "\"%s\"\n"
                   " does not correspond to any boundary face.\n"),
               missing, cpyfstr);
  }

  BFT_FREE(cpyfstr);
}

/*----------------------------------------------------------------------------
 * Give the list of the faces checking the constraint defined by
 * the Fortran string "fstr"
 *----------------------------------------------------------------------------*/

void CS_PROCF(csgfac, CSGFAC)
(
 const char          *const fstr,      /* <-- Fortran string */
 int                 *const len,       /* <-- String Length  */
 int                 *const n_faces,   /* --> number of faces */
 int                 *const faces      /* --> faces  */
 CS_ARGF_SUPP_CHAINE
)
{
  char * cpyfstr;
  int i, c_id;
  int lenwithoutblank;

  /* Initialization */
  *n_faces = 0;

  /* Copy fstr without last blanks  */
  lenwithoutblank = *len - 1;
  while(fstr[lenwithoutblank--] == ' ' &&  lenwithoutblank >= 0);

  if (lenwithoutblank < -1)
    return;

  lenwithoutblank += 2;

  BFT_MALLOC(cpyfstr, lenwithoutblank + 1, char);

  for(i = 0 ; i < lenwithoutblank; i++)
    cpyfstr[i] = fstr[i];
  cpyfstr[lenwithoutblank] = '\0';

  /* Get faces with C string */

  c_id = fvm_selector_get_list(cs_glob_mesh->select_i_faces,
                               cpyfstr,
                               n_faces,
                               faces);

  if (fvm_selector_n_missing(cs_glob_mesh->select_i_faces, c_id) > 0) {
    const char *missing
      = fvm_selector_get_missing(cs_glob_mesh->select_i_faces, c_id, 0);
    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_("The group or attribute \"%s\" in the\n"
                 "selection criterion:\n"
                 "\"%s\"\n"
                 " does not correspond to any interior face.\n"),
               missing, cpyfstr);
  }

  BFT_FREE(cpyfstr);
}

/*----------------------------------------------------------------------------
 * Give the list of the faces checking the constraint defined by
 * the Fortran string "fstr"
 *----------------------------------------------------------------------------*/

void CS_PROCF(csgcel, CSGCEL)
(
 const char          *const fstr,     /* <-- Fortran string */
 int                 *const len,      /* <-- String Length  */
 int                 *const n_cells,  /* --> number of cells */
 int                 *const cells     /* --> cells  */
 CS_ARGF_SUPP_CHAINE
)
{
  char * cpyfstr;
  int i, c_id;
  int lenwithoutblank;

  /* Initialization */
  *n_cells = 0;

  /* Copy fstr without last blanks  */
  lenwithoutblank = *len - 1;
  while(fstr[lenwithoutblank--] == ' ' &&  lenwithoutblank >= 0);

  if (lenwithoutblank < -1)
    return;

  lenwithoutblank += 2;

  BFT_MALLOC(cpyfstr, lenwithoutblank + 1, char);

  for(i = 0 ; i < lenwithoutblank; i++)
    cpyfstr[i] = fstr[i];
  cpyfstr[lenwithoutblank] = '\0';

  /* Get cells with C string */
  c_id = fvm_selector_get_list(cs_glob_mesh->select_cells,
                               cpyfstr,
                               n_cells,
                               cells);

  if (fvm_selector_n_missing(cs_glob_mesh->select_cells, c_id) > 0) {
    const char *missing
      = fvm_selector_get_missing(cs_glob_mesh->select_cells, c_id, 0);
    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_("The group or attribute \"%s\" in the selection\n"
                 "criteria:\n"
                 "\"%s\"\n"
                 " does not correspond to any cell.\n"),
               missing, cpyfstr);
  }

  BFT_FREE(cpyfstr);
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
