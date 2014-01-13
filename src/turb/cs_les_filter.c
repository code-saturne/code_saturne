/*============================================================================
 * Filters for dynamic models.
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <float.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_les_filter.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definition
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute filters for dynamic models. This function deals with the standard
 * or extended neighborhood.
 *
 * Fortran Interface :
 *
 * subroutine cfiltr (var, f_var, wbuf1, wbuf2)
 * *****************
 *
 * double precision(*) var[]   <-- array of variables to filter
 * double precision(*) f_var[] --> filtered variable array
 * double precision(*) wbuf1[] --- working buffer
 * double precision(*) wbuf2[] --- working buffer
 *----------------------------------------------------------------------------*/

void
CS_PROCF (cfiltr, CFILTR)(cs_real_t  var[],
                          cs_real_t  f_var[],
                          cs_real_t  wbuf1[],
                          cs_real_t  wbuf2[])
{
  cs_int_t  i, j, k;

  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_int_t  n_cells = mesh->n_cells;
  const cs_int_t  n_cells_ext = mesh->n_cells_with_ghosts;
  const cs_int_t  *cell_cells_idx = mesh->cell_cells_idx;
  const cs_int_t  *cell_cells_lst = mesh->cell_cells_lst;
  const cs_real_t  *cell_vol = cs_glob_mesh_quantities->cell_vol;

  assert(cell_cells_idx != NULL);

  /* Synchronize variable */

  if (mesh->halo != NULL)
    cs_halo_sync_var(mesh->halo, CS_HALO_EXTENDED, var);

  /* Allocate and initialize working buffers */

  for (i = 0; i < n_cells_ext; i++) {
    wbuf1[i] = 0;
    wbuf2[i] = 0;
  }

  /* Define filtered variable array */

  for (i = 0; i < n_cells; i++) {

    wbuf1[i] += var[i] * cell_vol[i];
    wbuf2[i] += cell_vol[i];

    /* Loop on connected cells (without cells sharing a face) */

    for (j = cell_cells_idx[i] - 1; j < cell_cells_idx[i+1] - 1; j++) {

      k = cell_cells_lst[j] - 1;
      wbuf1[i] += var[k] * cell_vol[k];
      wbuf2[i] += cell_vol[k];

    }

  } /* End of loop on cells */

  for (k = 0; k < mesh->n_i_faces; k++) {

    i = mesh->i_face_cells[2*k] - 1;
    j = mesh->i_face_cells[2*k + 1] - 1;

    wbuf1[i] += var[j] * cell_vol[j];
    wbuf2[i] += cell_vol[j];
    wbuf1[j] += var[i] * cell_vol[i];
    wbuf2[j] += cell_vol[i];

  }

  for (i = 0; i < n_cells; i++)
    f_var[i] = wbuf1[i]/wbuf2[i];

  /* Synchronize variable */

  if (mesh->halo != NULL)
    cs_halo_sync_var(mesh->halo, CS_HALO_STANDARD, f_var);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
