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

/*============================================================================
 * Functions dealing with the selection of cs_mesh_t entities
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <string.h>

#if defined(_CS_HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_error.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_mesh_select.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definition
 *============================================================================*/

struct _cs_mesh_select_t {

  cs_int_t     n_colors;    /* Number of colors assoicated to the selection  */
  cs_int_t    *colors;      /*  Color list */

  cs_int_t     n_groups;    /* Number of groups  */
  char       **groups;      /* Group name list */

  cs_bool_t    inv_selection;

};

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Creation and initialization of a cs_mesh_select_t structure.
 *
 * parameters:
 *   n_colors    --> number of colors
 *   n_groups    --> number of groups
 *   colors      --> color list
 *   groups      --> group name list
 *   invsel      --> invert selection if true
 *
 * returns:
 *   pointer to created cs_mesh_select_t structure
 *----------------------------------------------------------------------------*/

cs_mesh_select_t *
cs_mesh_select_create(cs_int_t     n_colors,
                      cs_int_t     n_groups,
                      cs_int_t    *colors,
                      char       **groups,
                      cs_bool_t    invsel)
{
  cs_int_t color_id, grp_id, length;

  cs_mesh_select_t *selection = NULL;

  BFT_MALLOC(selection, 1, cs_mesh_select_t);

  selection->inv_selection = invsel;

  /* Initialize colors */

  selection->n_colors = n_colors;
  BFT_MALLOC(selection->colors, n_colors, cs_int_t);

  for (color_id = 0; color_id < n_colors; color_id++)
    selection->colors[color_id] = colors[color_id];

  /* Initialize groups */

  selection->n_groups = n_groups;
  BFT_MALLOC(selection->groups, n_groups, char *);

  for (grp_id = 0; grp_id < n_groups; grp_id++) {

    length = strlen(groups[grp_id]);
    BFT_MALLOC(selection->groups[grp_id], length + 1, char);
    strncpy(selection->groups[grp_id], groups[grp_id], length);
    (selection->groups[grp_id])[length] = '\0';

  }

  return selection;
}

/*----------------------------------------------------------------------------
 * Destroy a cs_mesh_select_t structure
 *
 * parameters:
 *   selection --> pointer to selection structure that should be destroyed
 *
 * returns:
 *   NULL
 *----------------------------------------------------------------------------*/

cs_mesh_select_t *
cs_mesh_select_destroy(cs_mesh_select_t  *selection)
{
  cs_int_t grp_id;

  BFT_FREE(selection->colors);

  for (grp_id = 0; grp_id < selection->n_groups; grp_id++)
    BFT_FREE(selection->groups[grp_id]);
  BFT_FREE(selection->groups);

  BFT_FREE(selection);

  return NULL;
}

/*----------------------------------------------------------------------------
 * Get the number of colors in a cs_mesh_select_t structure.
 *
 * parameters:
 *   selection --> pointer to selection structure
 *
 * returns:
 *   Number of colors in the cs_mesh_select_t struture
 *----------------------------------------------------------------------------*/

cs_int_t
cs_mesh_select_get_n_colors(const cs_mesh_select_t *selection)
{

  if (selection == NULL)
    return 0;
  else
    return selection->n_colors;

}

/*----------------------------------------------------------------------------
 * Get the number of groups in a cs_mesh_select_t struture.
 *
 * parameters:
 *   selection --> pointer to a selection structure
 *
 * returns:
 *   Number of groups in the cs_mesh_select_t struture
 *----------------------------------------------------------------------------*/

cs_int_t
cs_mesh_select_get_n_groups(const cs_mesh_select_t *selection)
{

  if (selection == NULL)
    return 0;
  else
    return selection->n_groups;

}

/*----------------------------------------------------------------------------
 * Extract border faces from criteria included in a cs_mesh_select_t
 * structure.
 *
 * parameters:
 *   mesh               --> pointer to a mesh structure
 *   selection          --> pointer to a selection structure
 *   n_selected_b_faces <-- number of border faces selected
 *   b_face_select_lst  <-- list of the selected border faces
 *----------------------------------------------------------------------------*/

void
cs_mesh_select_extract_b_faces(const cs_mesh_t        *const mesh,
                               const cs_mesh_select_t       *selection,
                               cs_int_t                     *n_selected_b_faces,
                               cs_int_t                    **b_face_select_lst)
{
  cs_int_t  fam_id, id, idx, color_id, grp_id, i, grp_num ;
  cs_bool_t  end;

  char  *grp_name = NULL;
  cs_bool_t  *selected_families = NULL;
  cs_int_t   *_selected_lst = NULL;

  cs_int_t  *family_item = mesh->family_item;
  cs_int_t  *group_idx = mesh->group_idx;
  char      *group_lst = mesh->group_lst;

  BFT_MALLOC(selected_families, mesh->n_families, cs_bool_t);

  for (fam_id = 0; fam_id < mesh->n_families; fam_id++)
    selected_families[fam_id] = false;

  /* Look for families with criteria present in selection structure */

  for (id = 0, idx = 0; id < mesh->n_max_family_items; id++) {

    for (fam_id = 0; fam_id < mesh->n_families; fam_id++) {

      if (family_item[idx] > 0) { /* Color */

        color_id = 0;
        end = false;
        while (end == false && color_id < selection->n_colors) {

          if (family_item[idx] == selection->colors[color_id]) {
            selected_families[fam_id] = true;
            end = true;
          }
          color_id++;

        }

      } /* End of color treatment */

      else if (family_item[idx] < 0) { /* Group  */

        grp_id = 0;
        end = false;
        grp_num = CS_ABS(family_item[idx]);
        grp_name = &(group_lst[group_idx[grp_num - 1] - 1]);

        while (end == false && grp_id < selection->n_groups) {

          if (!strcmp(selection->groups[grp_id], grp_name)) {
            selected_families[fam_id] = true;
            end = true;
          }
          grp_id++;

        }

      } /* End of group treatment */

      else { /* Default family */

        assert(family_item[idx] == 0);

      }

      idx++;

    } /* End of loop on families */

  } /* End of loop on family properties */

  /* If inversion is required */

  if (selection->inv_selection == true) {
    for (fam_id = 0;  fam_id < mesh->n_families; fam_id++) {

      if (selected_families[fam_id] == false)
        selected_families[fam_id] = true;
      else if (selected_families[fam_id] == true)
        selected_families[fam_id] = false;

    }
  }

  /* Select border faces belonging to selected families */

  *n_selected_b_faces = 0;

  BFT_MALLOC(_selected_lst, mesh->n_b_faces, cs_int_t);

  for (i = 0; i < mesh->n_b_faces; i++) {
    if (selected_families[mesh->b_face_family[i] - 1] == true) {

      _selected_lst[*n_selected_b_faces] = i + 1;
      *n_selected_b_faces += 1;

    }
  }

  /* Realloc border face list if necessary */

  if (*n_selected_b_faces != mesh->n_b_faces)
    BFT_REALLOC(_selected_lst, *n_selected_b_faces, cs_int_t);

  *b_face_select_lst = _selected_lst;

  /* Free temporary buffers */

  BFT_FREE(selected_families);
}

/*----------------------------------------------------------------------------
 * Dump a cs_mesh_select_t structure.
 *
 * parameters:
 *   selection --> pointer to a selection structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_select_dump(cs_mesh_select_t  *selection)
{
  cs_int_t i;

  assert(selection != NULL);

  if (selection->inv_selection == true)
    bft_printf(_("Selection inversion requested\n"));

  bft_printf(_("\nColor(s) of coupled faces:\n"));
  bft_printf(_("Number of color(s): %i\n"),selection->n_colors);
  for (i = 0; i < selection->n_colors; i++)
    bft_printf(_("Color number: %i\n"), selection->colors[i]);

  bft_printf(_("\nGroup(s) of coupled faces:\n"));
  bft_printf(_("Number of group(s): %i\n"),selection->n_groups);
  for (i = 0; i < selection->n_groups; i++)
    bft_printf(_("Group: %s\n"), selection->groups[i]);

  bft_printf_flush();
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
