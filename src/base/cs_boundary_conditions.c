/*============================================================================
 * Boundary condition handling.
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

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_halo.h"
#include "cs_matrix.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_prototypes.h"
#include "cs_post.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_boundary_conditions.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/* Mappings to MPI datatypes */

#if defined(HAVE_MPI)

typedef struct
{
  int val;
  int rank;
} _mpi_int_int_t;

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Local Structure Definitions
 *----------------------------------------------------------------------------*/

/* Face marker structure for selecting error postprocessing output faces. */

typedef struct {

  cs_lnum_t       n_faces;   /* Number of boundary faces */
  unsigned char  *flag;      /* 0 for unmarked faces, 1 for marked faces */

} _error_face_marker_t;

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Transfer info on face with lowest global number to rank 0.
 *
 * parameters:
 *   face_gnum <-> pointer to global face number
 *   face_coo  <-> pointer to face coordinates
 *   face_type <-> pointer to face type
 *----------------------------------------------------------------------------*/

static void
_min_gnum_face(cs_gnum_t  *face_gnum,
               int        *face_type,
               double      face_coo[3])
{
#if defined(HAVE_MPI)

  /* local variables */

  cs_gnum_t  min_face_gnum;
  _mpi_int_int_t  val_in, val_min;

  /* Return immediately if not running under MPI */

  if (cs_glob_n_ranks < 2)
    return;

  /* Obtain the lowest global face number with an error; use minloc
     with a marker, rather than with a global number directly, in
     case cs_gnum_t is larger than an integer) */

  MPI_Allreduce(face_gnum, &min_face_gnum, 1, CS_MPI_GNUM, MPI_MIN,
                cs_glob_mpi_comm);

  if (*face_gnum == min_face_gnum)
    val_in.val = 0;
  else
    val_in.val = 1;
  val_in.rank = cs_glob_rank_id;

  MPI_Allreduce(&val_in, &val_min, 1, MPI_2INT, MPI_MINLOC,
                cs_glob_mpi_comm);

  /* Now exchange values */

  if (val_min.rank > 0) {

    if (val_min.rank == cs_glob_rank_id) {
      MPI_Send(face_gnum, 1, CS_MPI_GNUM, 0, 1, cs_glob_mpi_comm);
      MPI_Send(face_type, 1, MPI_INT, 0, 2, cs_glob_mpi_comm);
      MPI_Send(face_coo, 3, MPI_DOUBLE, 0, 3, cs_glob_mpi_comm);
    }
    else if (cs_glob_rank_id == 0) {
      MPI_Status status;
      MPI_Recv(face_gnum, 1, CS_MPI_GNUM, val_min.rank, 1,
               cs_glob_mpi_comm, &status);
      MPI_Recv(face_type, 1, MPI_INT, val_min.rank, 2,
               cs_glob_mpi_comm, &status);
      MPI_Recv(face_coo, 3, MPI_DOUBLE, val_min.rank, 3,
               cs_glob_mpi_comm, &status);
    }

  }

#endif
}

/*----------------------------------------------------------------------------
 * Function for selection of faces with boundary condition errors.
 *
 * parameters:
 *   input    <-- pointer to input (face marker structure)
 *   n_faces  --> number of selected faces
 *   face_ids --> array of selected face ids (0 to n-1 numbering)
 *----------------------------------------------------------------------------*/

static void
_post_error_faces_select(void         *input,
                         cs_lnum_t    *n_faces,
                         cs_lnum_t   **face_ids)
{
  cs_lnum_t face_id;

  cs_lnum_t _n_faces = 0;
  cs_lnum_t *_face_ids = NULL;

  const _error_face_marker_t  *marker = input;

  BFT_MALLOC(_face_ids, marker->n_faces, cs_lnum_t);

  for (face_id = 0; face_id < marker->n_faces; face_id++) {
    if (marker->flag[face_id] != 0)
      _face_ids[_n_faces++] = face_id;
  }

  *n_faces = _n_faces;
  *face_ids = _face_ids;
}

/*----------------------------------------------------------------------------
 * Function for selection of faces with valid boundary conditions.
 *
 * parameters:
 *   input    <-- pointer to input (face marker structure)
 *   n_faces  --> number of selected faces
 *   face_ids --> array of selected face ids (0 to n-1 numbering)
 *----------------------------------------------------------------------------*/

static void
_post_valid_faces_select(void         *input,
                         cs_lnum_t    *n_faces,
                         cs_lnum_t   **face_ids)
{
  cs_lnum_t face_id;

  cs_lnum_t _n_faces = 0;
  cs_lnum_t *_face_ids = NULL;

  const _error_face_marker_t  *marker = input;

  BFT_MALLOC(_face_ids, marker->n_faces, cs_lnum_t);

  for (face_id = 0; face_id < marker->n_faces; face_id++) {
    if (marker->flag[face_id] == 0)
      _face_ids[_n_faces++] = face_id;
  }

  *n_faces = _n_faces;
  *face_ids = _face_ids;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public Fortran function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Handling of boundary condition definition errors and associated output.
 *
 * For each boundary face, itypfb defines the boundary condition type.
 * As a convention here, zero values correspond to undefined types,
 * positive values to defined types (with no error), and negative values
 * to defined types with inconsistent or incompatible values, the
 * absolute value indicating the original boundary condition type.
 *
 * Fortran Interface:
 *
 * SUBROUTINE BCDERR (ITYPFB)
 * *****************
 *
 * INTEGER          ITYPFB      : <-> : Array of BC type ids
 *----------------------------------------------------------------------------*/

void CS_PROCF (bcderr, BCDERR)
(
 cs_int_t        *itypfb
)
{
  cs_boundary_conditions_error(itypfb);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Handling of boundary condition definition errors and associated output.
 *
 * For each boundary face, bc_type defines the boundary condition type.
 * As a convention here, zero values correspond to undefined types,
 * positive values to defined types (with no error), and negative values
 * to defined types with inconsistent or incompatible values, the
 * absolute value indicating the original boundary condition type.
 *
 *
 * parameters:
 *   bc_type   <-- array of BC type ids
 *----------------------------------------------------------------------------*/

void
cs_boundary_conditions_error(const cs_int_t  bc_type[])
{
  /* local variables */

  cs_lnum_t  face_id;
  _error_face_marker_t  marker;

  cs_gnum_t  n_errors = 0;

  const cs_mesh_t *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t *mesh_q = cs_glob_mesh_quantities;

  const cs_lnum_t n_b_faces = mesh->n_b_faces;

  /* Prepare face marker */

  marker.n_faces = n_b_faces;
  BFT_MALLOC(marker.flag, marker.n_faces, unsigned char);

  for (face_id = 0; face_id < n_b_faces; face_id++)
    marker.flag[face_id] = 0;

  /* Count and mark faces with problems */

  {
    int        err_face_type;
    cs_real_t  err_face_coo[3];
    cs_gnum_t  err_face_gnum = 0;

    const cs_int_t  *_bc_type = bc_type;

    for (face_id = 0; face_id < n_b_faces; face_id++) {

      /* _bc_type[] used to determine if we have an error */

      if (_bc_type[face_id] < 1) {

        cs_gnum_t face_gnum;

        if (mesh->global_b_face_num != NULL)
          face_gnum = mesh->global_b_face_num[face_id];
        else
          face_gnum = face_id + 1;

        marker.flag[face_id] = 1;

        if (err_face_gnum == 0 || face_gnum < err_face_gnum) {
          int coo_id;
          err_face_type = _bc_type[face_id];
          for (coo_id = 0; coo_id < 3; coo_id++)
            err_face_coo[coo_id] = mesh_q->b_face_cog[face_id*3 + coo_id];
        }

        n_errors += 1;
      }
    }

    cs_parall_counter(&n_errors, 1);

    /* Obtain the lowest global face number with an error,
       and print associated info */

    _min_gnum_face(&err_face_gnum, &err_face_type, err_face_coo);

    if (cs_glob_rank_id < 1) {

      bft_printf(_("\nFirst face with boundary condition definition error\n"
                   "  (out of %llu)\n"
                   "  has boundary condition type %d, center (%g, %g, %g)\n\n"),
                 (unsigned long long)n_errors , abs(err_face_type),
                 err_face_coo[0], err_face_coo[1], err_face_coo[2]);
    }

  }

  /* If post-processing is active, output boundary condition info */
  /*--------------------------------------------------------------*/

  if (mesh->b_face_vtx_idx) {

    int ii;

    cs_gnum_t n_valid_faces = 0;
    int mesh_id[2] = {0, 0};

    const int writer_id = -2;
    const int writer_ids[] = {writer_id};

    n_errors = 0;

    cs_post_init_error_writer();

    /* Mesh for invalid faces */

    mesh_id[0] = cs_post_get_free_mesh_id();

    cs_post_define_surface_mesh_by_func(mesh_id[0],
                                        _("Faces with B.C. error"),
                                        NULL,
                                        _post_error_faces_select,
                                        NULL,
                                        &marker,
                                        false, /* time varying */
                                        true,  /* add groups if present */
                                        false, /* auto variables */
                                        1,
                                        writer_ids);

    /* Mesh for valid faces */

    for (face_id = 0; face_id < n_b_faces; face_id++) {
      if (marker.flag[face_id] == 0)
        n_valid_faces += 1;
    }

    cs_parall_counter(&n_valid_faces, 1);

    if (n_valid_faces > 0) {

      mesh_id[1] = cs_post_get_free_mesh_id();

      cs_post_define_surface_mesh_by_func(mesh_id[1],
                                          _("Faces with valid B.C.'s"),
                                          NULL,
                                          _post_valid_faces_select,
                                          NULL,
                                          &marker,
                                          false, /* time varying */
                                          true,  /* add groups if present */
                                          false, /* auto variables */
                                          1,
                                          writer_ids);

    }

    cs_post_activate_writer(writer_id, 1);

    cs_post_write_meshes(NULL);

    BFT_FREE(marker.flag);

    {
      size_t name_size = 0;
      char var_name[32];

      const cs_int_t  *_bc_type = bc_type;

      var_name[0] = '\0';
      strncpy(var_name + name_size, _("BC type"), 31 - name_size);

      for (ii = 0; ii < 2; ii++) {

        if (mesh_id[ii] != 0)
          cs_post_write_var(mesh_id[ii],
                            var_name,
                            1,
                            false, /* no interlace */
                            true,  /* use parents */
                            CS_POST_TYPE_cs_int_t,
                            NULL,
                            NULL,
                            _bc_type,
                            NULL);

      }
    }
  }

  bft_error
    (__FILE__, __LINE__, 0,
     _("\nSome boundary condition definitions are incomplete or incorrect.\n\n"
       "  For details, read the end of the calculation log,\n"
       "  or visualize the error postprocessing output."));
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
