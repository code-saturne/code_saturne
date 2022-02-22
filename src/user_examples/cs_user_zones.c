/*============================================================================
 * Boundary and volume zone definitions.
 *============================================================================*/

/* VERS */

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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local (user defined) function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Example function for advanced selection of cells adjacent to a boundary
 * group named "G3".
 *
 * If non-empty and not containing all elements, a list of elements
 * of the parent mesh belonging to the location should be allocated
 * (using BFT_MALLOC) and defined by this function when called.
 * This list's lifecycle is then managed by the mesh location object.
 *
 * Note: if the input pointer is non-NULL, it must point to valid data
 * when the selection function is called, so that value or structure should
 * not be temporary (i.e. local);
 *
 * parameters:
 *   input       <-- pointer to optional (untyped) value or structure.
 *   m           <-- pointer to associated mesh structure.
 *   location_id <-- id of associated location.
 *   n_elts      --> number of selected elements
 *   elt_ids     --> list of selected elements (0 to n-1 numbering).
 *----------------------------------------------------------------------------*/

/*! [select_func_vol2] */
static void
_g3_boundary_cells(void              *input,
                   const cs_mesh_t   *m,
                   int                location_id,
                   cs_lnum_t         *n_elts,
                   cs_lnum_t        **elt_ids)
{
  CS_UNUSED(input);
  CS_UNUSED(location_id);

  /* Allocate selection list */

  cs_lnum_t n_b_faces = 0;
  cs_lnum_t *b_face_ids = NULL;

  BFT_MALLOC(b_face_ids, m->n_b_faces, cs_lnum_t);

  cs_selector_get_b_face_list("G3", &n_b_faces, b_face_ids);

  char *cell_flag;
  BFT_MALLOC(cell_flag, m->n_cells, char);

  for (cs_lnum_t i = 0; i < m->n_cells; i++)
    cell_flag[i] = 0;

  /* Now that flag is built, test for adjacency */

  for (cs_lnum_t i= 0; i < n_b_faces; i++) {
    cs_lnum_t face_id = b_face_ids[i];
    cs_lnum_t c_id = m->b_face_cells[face_id];
    cell_flag[c_id] = 1;
  }

  /* Now build list from flag */

  cs_lnum_t _n_elts = 0;
  for (cs_lnum_t i = 0; i < m->n_cells; i++) {
    if (cell_flag[i] != 0)
      _n_elts++;
  }

  cs_lnum_t *_elt_ids;
  BFT_MALLOC(_elt_ids, _n_elts, cs_lnum_t);

  _n_elts = 0;
  for (cs_lnum_t i = 0; i < m->n_cells; i++) {
    if (cell_flag[i] != 0) {
      _elt_ids[_n_elts] = i;
      _n_elts++;
    }
  }

  /* Free memory */

  BFT_FREE(b_face_ids);
  BFT_FREE(cell_flag);

  /* Set return values */

  *n_elts = _n_elts;
  *elt_ids = _elt_ids;
}
/*! [select_func_vol2] */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define volume and surface zones.
 *
 * See \ref sec_selection_criteria for details on selection criteria.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_zones(void)
{
  /* Example:
     a volume zone for head losses based on geometrical criteria
     described by a rule. */

  /*! [user_zones_head_loss_1] */
  {
    cs_volume_zone_define("head_loss_1",
                          "box[4.0, 6.0, -1e6, 2.0, 8.0, 1e6]",
                          CS_VOLUME_ZONE_HEAD_LOSS);
  }
  /*! [user_zones_head_loss_1] */

  /* Example:
     define 3 spherical source term zones based on geometrical criteria
     described by a rule. */

  /*! [user_zones_volume_1] */
  {
    char name[128], criteria[128];

    for (int i = 0; i < 3; i++) {

      double s_coords[] = {0, 0, 2.0*i};

      snprintf(name, 127, "source_%d", i);
      snprintf(criteria, 127, "sphere[%f, %f, %f, 0.5]",
               s_coords[0], s_coords[1], s_coords[2]);

      cs_volume_zone_define(name, criteria, CS_VOLUME_ZONE_SOURCE_TERM);
    }
  }
  /*! [user_zones_volume_1] */

  /*! [user_zones_volume_2] */
  {
    char name[128], criteria[128];

    for (int i = 0; i < 3; i++) {

      double s_coords[] = {0, 0, 2.0*i};

      snprintf(name, 127, "source_%d", i);
      snprintf(criteria, 127, "sphere[%f, %f, %f, 0.5]",
               s_coords[0], s_coords[1], s_coords[2]);

      cs_volume_zone_define_by_func("G3_B", _g3_boundary_cells, NULL, 0);
    }
  }
  /*! [user_zones_volume_2] */

   /* Example:
     Define zones from a text file defining turbines modelled
     with the actuator disk approach. */

  {
    char name[128], criteria[128];

    FILE* file = fopen("turbines", "rt");
    int n_turbines = 0;

    /* Some parameters */
    cs_real_t turb_lenght = 1.;
    cs_real_t radius = 126./2.;

    cs_real_t wind_dir = 0.;

    if (fscanf(file, "%d\n", &n_turbines) != 1)
      bft_error(__FILE__,__LINE__, 0,
                _("Could not read the number of turbines."));

    for (int i = 0; i < n_turbines; i++) {
      int turbine_id = 0;
      if (fscanf(file, "%d\n", &turbine_id) != 1)
        bft_error(__FILE__,__LINE__, 0, _("Could not read turbine %d."), i);

      float turb_cen[3];
      if (fscanf(file, "%f\n", &(turb_cen[0])) != 1)
        bft_error(__FILE__,__LINE__, 0, _("Could not read turbine x."));
      if (fscanf(file, "%f\n", &(turb_cen[1])) != 1)
        bft_error(__FILE__,__LINE__, 0, _("Could not read turbine y."));
      if (fscanf(file, "%f\n", &(turb_cen[2])) != 1)
        bft_error(__FILE__,__LINE__, 0, _("Could not read turbine z."));

      double s_coords[]
        = {turb_cen[0] - 0.5 * turb_lenght * cos(wind_dir),
           turb_cen[1] - 0.5 * turb_lenght * sin(wind_dir),
           turb_cen[2]};
      double e_coords[]
        = {turb_cen[0] + 0.5 * turb_lenght * cos(wind_dir),
           turb_cen[1] + 0.5 * turb_lenght * sin(wind_dir),
           turb_cen[2]};

      snprintf(name, 127, "turbine_%d", i);
      snprintf(criteria, 127, "cylinder[%f, %f, %f, %f, %f, %f, %f]",
               s_coords[0], s_coords[1], s_coords[2], /* Start */
               e_coords[0], e_coords[1], e_coords[2], /* End */
               radius); /* Radius */

      cs_volume_zone_define(name, criteria, CS_VOLUME_ZONE_SOURCE_TERM);
    }
    if (fclose(file) != 0)
      bft_error(__FILE__,__LINE__, 0, _("Could not close the file."));
  }

  /* Example:
     define simple boundary zones, allowing all faces not in the
     "INLET" or "OUTLET" groups to be considered as walls */

  /*! [user_zones_boundary_1] */
  {
    int z_id = cs_boundary_zone_define("wall", "all[]", 0);
    cs_boundary_zone_set_overlay(z_id, true);

    cs_boundary_zone_define("inlet", "INLET", 0);
    cs_boundary_zone_define("outlet", "OUTLET", 0);
  }
  /*! [user_zones_boundary_1] */

  /* Example:
     allow overlapping of an existing zone named "wall"
   */

  /*! [user_zones_boundary_2] */
  {
    int z_id = cs_boundary_zone_by_name("wall")->id;
    cs_boundary_zone_set_overlay(z_id, true);
  }
  /*! [user_zones_boundary_2] */

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
