/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2009 EDF S.A., France
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
 * Boundary condition handling.
 *============================================================================*/

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
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_error.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <fvm_defs.h>
#include <fvm_parall.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_halo.h"
#include "cs_matrix.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_prototypes.h"
#include "cs_post.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_boundary_conditions.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

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
_min_gnum_face(fvm_gnum_t  *face_gnum,
               int         *face_type,
               double       face_coo[3])
{
#if defined(HAVE_MPI)

  /* local variables */

  fvm_gnum_t  min_face_gnum;
  _mpi_int_int_t  val_in, val_min;

  /* Return immediately if not running under MPI */

  if (cs_glob_n_ranks < 2)
    return;

  /* Obtain the lowest global face number with an error; use minloc
     with a marker, rather than with a global number directly, in
     case fvm_gnum_t is larger than an integer) */

  MPI_Allreduce(face_gnum, &min_face_gnum, 1, FVM_MPI_GNUM, MPI_MIN,
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
      MPI_Send(face_gnum, 1, FVM_MPI_GNUM, 0, 1, cs_glob_mpi_comm);
      MPI_Send(face_type, 1, MPI_INT, 0, 2, cs_glob_mpi_comm);
      MPI_Send(face_coo, 3, MPI_DOUBLE, 0, 3, cs_glob_mpi_comm);
    }
    else if (cs_glob_rank_id == 0) {
      MPI_Status status;
      MPI_Recv(face_gnum, 1, FVM_MPI_GNUM, val_min.rank, 1,
               cs_glob_mpi_comm, &status);
      MPI_Recv(face_type, 1, MPI_INT, val_min.rank, 2,
               cs_glob_mpi_comm, &status);
      MPI_Recv(face_coo, 3, MPI_DOUBLE, val_min.rank, 3,
               cs_glob_mpi_comm, &status);
    }

  }

#endif
}

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

  fvm_lnum_t  face_id;

  fvm_gnum_t  n_errors = 0;

  unsigned char  *face_marker = NULL;

  const cs_mesh_t *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t *mesh_q = cs_glob_mesh_quantities;

  const fvm_lnum_t n_b_faces = mesh->n_b_faces;

  /* Prepare face marker */

  BFT_MALLOC(face_marker, n_b_faces, unsigned char);

  for (face_id = 0; face_id < n_b_faces; face_id++)
    face_marker[face_id] = 0;

  /* Count and mark faces with problems */

  {
    int         err_face_type;
    cs_real_t   err_face_coo[3];
    fvm_gnum_t  err_face_gnum = 0;

    const cs_int_t  *_bc_type = bc_type;

    for (face_id = 0; face_id < n_b_faces; face_id++) {

      /* _bc_type[] used to determine if we have an error */

      if (_bc_type[face_id] < 1) {

        fvm_gnum_t face_gnum;

        if (mesh->global_b_face_num != NULL)
          face_gnum = mesh->global_b_face_num[face_id];
        else
          face_gnum = face_id + 1;

        face_marker[face_id] = 1;

        if (err_face_gnum == 0 || face_gnum < err_face_gnum) {
          int coo_id;
          err_face_type = _bc_type[face_id];
          for (coo_id = 0; coo_id < 3; coo_id++)
            err_face_coo[coo_id] = mesh_q->b_face_cog[face_id*3 + coo_id];
        }

        n_errors += 1;
      }
    }

    fvm_parall_counter(&n_errors, 1);

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

  if (mesh->i_face_vtx_idx != NULL || mesh->b_face_vtx_idx) {

    int ii;

    cs_int_t face_list_size = 0;
    cs_int_t *face_list = NULL;
    int mesh_id[2] = {0, 0};

    const int writer_id = -2;
    const int writer_ids[] = {writer_id};

    n_errors = 0;

    cs_post_init_error_writer();

    /* Prepare face marker */

    BFT_MALLOC(face_list, n_b_faces, cs_int_t);

    /* Mesh for invalid faces */

    face_list_size = 0;
    for (face_id = 0; face_id < n_b_faces; face_id++) {
      if (face_marker[face_id] != 0)
        face_list[face_list_size++] = face_id + 1;
    }

    mesh_id[0] = cs_post_get_free_mesh_id();

    cs_post_define_surface_mesh_by_list(mesh_id[0],
                                        _("Faces with B.C. error"),
                                        0,
                                        face_list_size,
                                        NULL,
                                        face_list,
                                        true,  /* add groups if present */
                                        false, /* auto variables */
                                        1,
                                        writer_ids);

    /* Mesh for valid faces */

    face_list_size = 0;
    for (face_id = 0; face_id < n_b_faces; face_id++) {
      if (face_marker[face_id] == 0)
        face_list[face_list_size++] = face_id + 1;
    }

    n_errors = face_list_size;
    fvm_parall_counter(&n_errors, 1);

    if (n_errors < mesh->n_g_b_faces) {

      mesh_id[1] = cs_post_get_free_mesh_id();

      cs_post_define_surface_mesh_by_list(mesh_id[1],
                                          _("Faces with valid B.C.'s"),
                                          0,
                                          face_list_size,
                                          NULL,
                                          face_list,
                                          true,  /* add groups if present */
                                          false, /* auto variables */
                                          1,
                                          writer_ids);

    }

    BFT_FREE(face_marker);

    cs_post_activate_writer(writer_id, 1);

    cs_post_write_meshes(-1, 0.0);

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
                            -1,
                            0.0,
                            NULL,
                            NULL,
                            _bc_type);

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
