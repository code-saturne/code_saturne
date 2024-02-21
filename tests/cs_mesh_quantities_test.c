/*============================================================================
 * Unit test for some mesh quantities algorithms.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <bft_error.h>
#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_math.h"
#include "cs_mesh_quantities.h"

/*---------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Face_quanties computation for non-convex face.
 ----------------------------------------------------------------------------*/

static void
_non_convex_face_quantities(void)
{
  /*
    We test the computation of face quantities on the following
    face, for which the center of gravity of vertices is outside
    the face.

    {0., 0.8}    {0.2, 0.8}
         --------
         |      |
         |      |
         |      |
         |      |------------------------------------  {1, 0.2}
         |                                          |
         --------------------------------------------
    {0., 0}                                            {1., 0.}

  */

  cs_real_3_t vtx_coord[]
    = {{0, 0, 0},
       {1, 0, 0},
       {1, 0.2, 0},
       {0.2, 0.2, 0},
       {0.2, 0.8, 0},
       {0., 0.8, 0}};
  cs_lnum_t face_vtx_idx[] = {0, 6};
  cs_lnum_t face_vtx_lst[] = {0, 1, 2, 3, 4, 5};

  cs_real_t face_cog[3], face_normal[3];

  cs_mesh_quantities_compute_face_quantities
    (1,
     vtx_coord,
     face_vtx_idx,
     face_vtx_lst,
     &face_cog,
     &face_normal);

  double w_a = 1.*0.2;
  double w_b = 0.2*(0.8-0.2);
  double w_tot = w_a + w_b;

  printf("%s:\n"
         "  expected COG: {%g, %g, %g}; surface: %g\n"
         "  computed COG: {%g %g, %g}; surface: %g\n",
         __func__,
         ((0+1)*0.5*w_a + (0+0.2)*0.5*w_b) / w_tot,
         ((0+0.2)*0.5*w_a + (0.2+0.8)*0.5*w_b) / w_tot,
         0.0,
         w_a + w_b,
         face_cog[0], face_cog[1], face_cog[2], face_normal[2]);
}

/*---------------------------------------------------------------------------*/

int
main (int argc, char *argv[])
{
  CS_UNUSED(argc);
  CS_UNUSED(argv);

  _non_convex_face_quantities();

  exit(EXIT_SUCCESS);
}
