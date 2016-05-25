/*============================================================================
 * Boundary condition handling.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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

#include <ple_locator.h>

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_coupling.h"
#include "cs_gradient.h"
#include "cs_field.h"
#include "cs_field_operator.h"
#include "cs_halo.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_connect.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_prototypes.h"
#include "cs_post.h"
#include "fvm_nodal.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_boundary_conditions.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_boundary_conditions.c
        Boundary condition handling.
*/

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

static int *_bc_type;

const int *cs_glob_bc_type = NULL;

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_boundary_conditions_mapped_set(int                        field_id,
                                    ple_locator_t             *locator,
                                    cs_mesh_location_type_t    location_type,
                                    int                        enforce_balance,
                                    int                        interpolate,
                                    cs_lnum_t                  n_faces,
                                    const cs_lnum_t           *faces,
                                    cs_real_t                 *balance_w,
                                    int                        nvarcl,
                                    cs_real_t                 *rcodcl);

void
cs_f_boundary_conditions_type_get_pointer(int  **itypfb);

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

/*----------------------------------------------------------------------------
 * Compute balance at inlet
 *
 * parameters:
 *   var_id          <-- variable id
 *   f               <-- associated field
 *   m               <-- associated mesh
 *   mq              <-- mesh quantities
 *   enforce_balance <-- balance handling option:
 *                         0: values are simply mapped
 *                         1: values are mapped, then multiplied
 *                            by a constant factor so that their
 *                            surface integral on selected faces
 *                            is preserved (relative to the
 *                            input values)
 *                         2: as 1, but with a boundary-defined
 *                            weight, defined by balance_w
 *                         3: as 1, but with a cell-defined
 *                            weight, defined by balance_w
 *   n_faces         <-- number of selected boundary faces
 *   faces           <-- list of selected boundary faces (0 to n-1),
 *                       or NULL if no indirection is needed
 *   balance_w       <-- optional balance weight, or NULL
 *   nvarcl          <-- number of variables requiring BC's
 *   rcodcl          <-- boundary condition values
 *   inlet_sum       --> inlet sum
 *----------------------------------------------------------------------------*/

static void
_inlet_sum(int                          var_id,
           const cs_field_t            *f,
           const cs_mesh_t             *m,
           const cs_mesh_quantities_t  *mq,
           int                          enforce_balance,
           cs_lnum_t                    n_faces,
           const cs_lnum_t             *faces,
           cs_real_t                   *balance_w,
           int                          nvarcl,
           cs_real_t                    rcodcl[],
           cs_real_t                    inlet_sum[])
{
  const int dim = f->dim;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_real_t *f_surf = mq->b_f_face_surf;

  /* Get field's variable id */

  assert(dim <= 9);

  for (cs_lnum_t j = 0; j < dim; j++) {

    cs_real_t  *_rcodcl = rcodcl + (var_id+j)*n_b_faces;

    inlet_sum[j] = 0.;

    if (enforce_balance == 1) {
      for (cs_lnum_t i = 0; i < n_faces; i++) {
        const cs_lnum_t f_id = (faces != NULL) ? faces[i] : i;
        inlet_sum[j] +=_rcodcl[f_id]*f_surf[f_id];
      }
    }
    else if (enforce_balance == 2) {
      for (cs_lnum_t i = 0; i < n_faces; i++) {
        const cs_lnum_t f_id = (faces != NULL) ? faces[i] : i;
        inlet_sum[j] +=_rcodcl[f_id]*f_surf[f_id]*balance_w[f_id];
      }
    }
    else if (enforce_balance == 3) {
      for (cs_lnum_t i = 0; i < n_faces; i++) {
        const cs_lnum_t f_id = (faces != NULL) ? faces[i] : i;
        const cs_lnum_t c_id = m->b_face_cells[f_id];
        inlet_sum[j] +=_rcodcl[f_id]*f_surf[f_id]*balance_w[c_id];
      }
    }

  }

  cs_parall_sum(dim, CS_REAL_TYPE, inlet_sum);
}

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Set mapped boundary conditions for a given field and mapping locator.
 *
 * parameters:
 *   field_id        <-- id of field whose boundary conditions are set
 *   locator         <-- associated mapping locator, as returned
 *                       by cs_boundary_conditions_map().
 *   location_type   <-- matching values location (CS_MESH_LOCATION_CELLS or
 *                       CS_MESH_LOCATION_BOUNDARY_FACES)
 *   enforce_balance <-- balance handling option:
 *                         0: values are simply mapped
 *                         1: values are mapped, then multiplied
 *                            by a constant factor so that their
 *                            surface integral on selected faces
 *                            is preserved (relative to the
 *                            input values)
 *                         2: as 1, but with a boundary-defined
 *                            weight, defined by balance_w
 *                         3: as 1, but with a cell-defined
 *                            weight, defined by balance_w
 *   interpolate     <-- interpolation option:
 *                         0: values are simply based on matching
 *                            cell or face center values
 *                         1: values are based on matching cell
 *                            or face center values, corrected
 *                            by gradient interpolation
 *   n_faces         <-- number of selected boundary faces
 *   faces           <-- list of selected boundary faces (1 to n),
 *                       or NULL if no indirection is needed
 *   balance_w       <-- optional balance weight, or NULL
 *   nvarcl          <-- number of variables requiring BC's
 *   rcodcl          <-> boundary condition values
 *----------------------------------------------------------------------------*/

void
cs_f_boundary_conditions_mapped_set(int                        field_id,
                                    ple_locator_t             *locator,
                                    cs_mesh_location_type_t    location_type,
                                    int                        enforce_balance,
                                    int                        interpolate,
                                    cs_lnum_t                  n_faces,
                                    const cs_lnum_t           *faces,
                                    cs_real_t                 *balance_w,
                                    int                        nvarcl,
                                    cs_real_t                 *rcodcl)
{
  cs_lnum_t *_faces = NULL;

  if (faces != NULL) {
    BFT_MALLOC(_faces, n_faces, cs_lnum_t);
    for (cs_lnum_t i = 0; i < n_faces; i++)
      _faces[i] = faces[i] - 1;
  }

  cs_boundary_conditions_mapped_set(cs_field_by_id(field_id),
                                    locator,
                                    location_type,
                                    enforce_balance,
                                    interpolate,
                                    n_faces,
                                    _faces,
                                    balance_w,
                                    nvarcl,
                                    rcodcl);

  BFT_FREE(_faces);
}

/*----------------------------------------------------------------------------
 * Get pointer to cs_glob_bc_type
 *----------------------------------------------------------------------------*/

void
cs_f_boundary_conditions_type_get_pointer(int **itypfb)
{
  *itypfb = _bc_type;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Handling of boundary condition definition errors and
 *        associated output.
 *
 * This function checks for errors, and simply returns if no errors are
 * encountered. In case of error, it outputs helpful information so as to
 * make it easier to locate the matching faces.
 *
 * For each boundary face, bc_type defines the boundary condition type.
 * As a convention here, zero values correspond to undefined types,
 * positive values to defined types (with no error), and negative values
 * to defined types with inconsistent or incompatible values, the
 * absolute value indicating the original boundary condition type.
 *
 * An optional label may be used if the error is related to another
 * attribute than the boundary type, for appropriate error reporting.
 *
 * \param[in]  bc_type     array of BC type ids
 * \param[in]  type_name   name of attribute in error, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_error(const int       *bc_flag,
                             const char      *type_name)
{
  /* local variables */

  cs_lnum_t  face_id;
  _error_face_marker_t  marker;

  cs_gnum_t  n_errors = 0;

  const cs_mesh_t *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t *mesh_q = cs_glob_mesh_quantities;

  const cs_lnum_t n_b_faces = mesh->n_b_faces;

  const char type_name_default[] = N_("boundary condition type");
  const char *_type_name = type_name_default;

  if (type_name != NULL)
    _type_name = type_name;

  /* First, simply check for faces with problems;
     _bc_flag[] used to determine if we have an error */

  const cs_int_t  *_bc_flag = bc_flag;

  for (face_id = 0; face_id < n_b_faces; face_id++) {
    if (_bc_flag[face_id] < 1)
      n_errors += 1;
  }

  cs_parall_counter(&n_errors, 1);

  if (n_errors < 1)
    return;

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

    for (face_id = 0; face_id < n_b_faces; face_id++) {

      /* _bc_flag[] used to determine if we have an error */

      if (_bc_flag[face_id] < 1) {

        cs_gnum_t face_gnum;

        if (mesh->global_b_face_num != NULL)
          face_gnum = mesh->global_b_face_num[face_id];
        else
          face_gnum = face_id + 1;

        marker.flag[face_id] = 1;

        if (err_face_gnum == 0 || face_gnum < err_face_gnum) {
          int coo_id;
          err_face_type = _bc_flag[face_id];
          for (coo_id = 0; coo_id < 3; coo_id++)
            err_face_coo[coo_id] = mesh_q->b_face_cog[face_id*3 + coo_id];
        }

      }
    }

    /* Obtain the lowest global face number with an error,
       and print associated info */

    _min_gnum_face(&err_face_gnum, &err_face_type, err_face_coo);

    if (cs_glob_rank_id < 1) {

      bft_printf(_("\nFirst face with boundary condition definition error\n"
                   "  (out of %llu)\n"
                   "  has %s %d, center (%g, %g, %g)\n\n"),
                 (unsigned long long)n_errors,
                 _type_name, abs(err_face_type),
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
                            _bc_flag,
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
/*!
 * \brief Locate shifted boundary face coordinates on possibly filtered
 *        cells or boundary faces for later interpolation.
 *
 * \param[in]  location_type    matching values location (CS_MESH_LOCATION_CELLS
 *                              or CS_MESH_LOCATION_BOUNDARY_FACES)
 * \param[in]  n_location_elts  number of selected location elements
 * \param[in]  n_faces          number of selected boundary faces
 * \param[in]  location_elts    list of selected location elements (0 to n-1),
 *                              or NULL if no indirection is needed
 * \param[in]  faces            list of selected boundary faces (0 to n-1),
 *                              or NULL if no indirection is needed
 * \param[in]  coord_shift      array of coordinates shift relative to selected
 *                              boundary faces
 * \param[in]  coord_stride     access stride in coord_shift: 0 for uniform
 *                              shift, 1 for "per face" shift.
 * \param[in]  tolerance        relative tolerance for point location.
 *
 * \return  associated locator structure
 */
/*----------------------------------------------------------------------------*/

ple_locator_t *
cs_boundary_conditions_map(cs_mesh_location_type_t    location_type,
                           cs_lnum_t                  n_location_elts,
                           cs_lnum_t                  n_faces,
                           const cs_lnum_t           *location_elts,
                           const cs_lnum_t           *faces,
                           cs_real_3_t               *coord_shift,
                           int                        coord_stride,
                           double                     tolerance)
{
  ple_locator_t *locator = NULL;

  /* Build temporary "donor" location  mesh */
  /*----------------------------------------*/

  fvm_nodal_t *nm = NULL;

  /* Convert from 0-based to 1-based numbering (note that in practice,
     we probably use 1-based numbering upstream, and convert back to
     1-based downstream, but we prefer to use a 0-based API here so as
     to make progress on global 1-based to 0-based conversion). */

  cs_lnum_t *_location_elts = NULL;
  if (location_elts != NULL) {
    BFT_MALLOC(_location_elts, n_location_elts, cs_lnum_t);
    for (cs_lnum_t i = 0; i < n_location_elts; i++)
      _location_elts[i] = location_elts[i] + 1;
  }

  if (location_type == CS_MESH_LOCATION_CELLS)
    nm = cs_mesh_connect_cells_to_nodal(cs_glob_mesh,
                                        "search mesh",
                                        false, /* include_families */
                                        n_location_elts,
                                        _location_elts);
  else if (location_type == CS_MESH_LOCATION_BOUNDARY_FACES)
    nm = cs_mesh_connect_faces_to_nodal(cs_glob_mesh,
                                        "search mesh",
                                        false, /* include_families */
                                        0,
                                        n_location_elts,
                                        NULL,
                                        _location_elts);

  BFT_FREE(_location_elts);

  /* Now build locator */
  /*-------------------*/

#if defined(PLE_HAVE_MPI)
  locator = ple_locator_create(cs_glob_mpi_comm, cs_glob_n_ranks, 0);
#else
  locator = ple_locator_create();
#endif

  int options[PLE_LOCATOR_N_OPTIONS];
  for (int i = 0; i < PLE_LOCATOR_N_OPTIONS; i++)
    options[i] = 0;
  options[PLE_LOCATOR_NUMBERING] = 0; /* base 0 numbering */

  /* Build location coordinates */

  ple_coord_t *point_coords;

  const cs_real_3_t *restrict b_face_cog
    = (const cs_real_3_t *restrict)cs_glob_mesh_quantities->b_face_cog;

  BFT_MALLOC(point_coords, n_faces*3, ple_coord_t);
  if (faces != NULL) {
    for (cs_lnum_t i = 0; i < n_faces; i++) {
      const cs_lnum_t face_id = faces[i];
      for (cs_lnum_t j = 0; j < 3; j++)
        point_coords[i*3 + j] =   b_face_cog[face_id][j]
                                + coord_shift[i*coord_stride][j];
    }
  }

  ple_locator_set_mesh(locator,
                       nm,
                       options,
                       0., /* tolerance_base */
                       tolerance,
                       3, /* dim */
                       n_faces,
                       NULL,
                       NULL, /* point_tag */
                       point_coords,
                       NULL, /* distance */
                       cs_coupling_mesh_extents,
                       cs_coupling_point_in_mesh_p);

  /* Check that location is OK */

  cs_gnum_t loc_count[2];
  loc_count[0] = ple_locator_get_n_exterior(locator);
  loc_count[1] = ple_locator_get_n_exterior(locator);
  cs_parall_counter(loc_count, 2);

  if (loc_count[1] > 0) {
    bft_error
      (__FILE__, __LINE__, 0,
       _("\nIn function %s,\n"
         "  %llu boundary faces (of %llu selected) were not matched to mesh\n"
         "  elements. Check your coordinate shift definitions."),
       __func__,
       (unsigned long long)loc_count[1],
       (unsigned long long)loc_count[0]);
  }

  BFT_FREE(point_coords);

  /* Shift from 1-base to 0-based locations */

  ple_locator_shift_locations(locator, -1);

  /* Nodal mesh is not needed anymore */

  nm = fvm_nodal_destroy(nm);

  /* Return initialized locator */

  return locator;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set mapped boundary conditions for a given field and mapping locator.
 *
 * \param[in]       f                field whose boundary conditions are set
 * \param[in]       locator          associated mapping locator, as returned
 *                                   by \ref cs_boundary_conditions_map.
 * \param[in]       location_type    matching values location
 *                                   (CS_MESH_LOCATION_CELLS or
 *                                   CS_MESH_LOCATION_BOUNDARY_FACES)
 * \param[in]       normalize        normalization option:
 *                                     0: values are simply mapped
 *                                     1: values are mapped, then multiplied
 *                                        by a constant factor so that their
 *                                        surface integral on selected faces
 *                                        is preserved (relative to the
 *                                        input values)
 *                                     2: as 1, but with a boundary-defined
 *                                        weight, defined by balance_w
 *                                     3: as 1, but with a cell-defined
 *                                       weight, defined by balance_w
 * \param[in]       interpolate      interpolation option:
 *                                     0: values are simply based on matching
 *                                        cell or face center values
 *                                     1: values are based on matching cell
 *                                        or face center values, corrected
 *                                        by gradient interpolation
 * \param[in]       n_faces          number of selected boundary faces
 * \param[in]       faces            list of selected boundary faces (0 to n-1),
 *                                   or NULL if no indirection is needed
 * \param[in]       balance_w        optional balance weight, or NULL
 * \param[in]       nvarcl           number of variables requiring BC's
 * \param[in, out]  rcodcl           boundary condition values
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_mapped_set(cs_field_t                *f,
                                  ple_locator_t             *locator,
                                  cs_mesh_location_type_t    location_type,
                                  int                        normalize,
                                  int                        interpolate,
                                  cs_lnum_t                  n_faces,
                                  const cs_lnum_t           *faces,
                                  cs_real_t                 *balance_w,
                                  int                        nvarcl,
                                  cs_real_t                  rcodcl[])
{
  int var_id = -1;

  const int dim = f->dim;
  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;
  const ple_lnum_t n_dist = ple_locator_get_n_dist_points(locator);
  const ple_lnum_t *dist_loc = ple_locator_get_dist_locations(locator);
  const ple_coord_t *dist_coords = ple_locator_get_dist_coords(locator);

  cs_field_interpolate_t   interpolation_type = CS_FIELD_INTERPOLATE_MEAN;

  cs_real_t inlet_sum_0[9], inlet_sum_1[9];
  cs_real_t *distant_var, *local_var;

  /* Get field's variable id */

  static int var_id_key = -1;
  if (var_id_key < 0)
    var_id_key = cs_field_key_id("variable_id");
  assert(var_id_key >= 0);

  var_id = cs_field_get_key_int(f, var_id_key) - 1;

  if (var_id < 0)
    return;

  assert(f->location_id == CS_MESH_LOCATION_CELLS);

  /* Compute initial balance if applicable */

  if (normalize > 0) {
    assert(dim <= 9);
    _inlet_sum(var_id,
               f,
               cs_glob_mesh,
               cs_glob_mesh_quantities,
               normalize,
               n_faces,
               faces,
               balance_w,
               nvarcl,
               rcodcl,
               inlet_sum_0);
  }

  /* Allocate working array */

  BFT_MALLOC(distant_var, n_dist*dim, cs_real_t);
  BFT_MALLOC(local_var, n_faces*dim, cs_real_t);

  /* Prepare values to send */
  /*------------------------*/

  if (interpolate)
    interpolation_type = CS_FIELD_INTERPOLATE_GRADIENT;

  assert(sizeof(ple_coord_t) == sizeof(cs_real_t));

  if (location_type == CS_MESH_LOCATION_CELLS || interpolate)
    cs_field_interpolate(f,
                         interpolation_type,
                         n_dist,
                         dist_loc,
                         (const cs_real_3_t *)dist_coords,
                         distant_var);

  else if (location_type == CS_MESH_LOCATION_BOUNDARY_FACES) {

    const cs_lnum_t *b_face_cells = cs_glob_mesh->b_face_cells;
    const cs_field_bc_coeffs_t   *bc_coeffs = f->bc_coeffs;

    /* If no boundary condition coefficients are available */

    if (bc_coeffs != NULL) {

      if (dim == 1) {
        for (cs_lnum_t i = 0; i < n_dist; i++) {
          cs_lnum_t f_id = dist_loc[i];
          cs_lnum_t c_id = b_face_cells[f_id];
          distant_var[i] =   bc_coeffs->a[f_id]
                           + bc_coeffs->b[f_id]*f->val[c_id];
        }
      }
      else {
        for (cs_lnum_t i = 0; i < n_dist; i++) {
          cs_lnum_t f_id = dist_loc[i];
          cs_lnum_t c_id = b_face_cells[f_id];
          for (cs_lnum_t j = 0; j < dim; j++) {
            distant_var[i*dim+j] = bc_coeffs->a[f_id*dim+j];
            for (cs_lnum_t k = 0; k < dim; k++)
              distant_var[i*dim+j] +=  bc_coeffs->b[(f_id*dim+k)*dim + j]
                                      *f->val[c_id*dim+k];
          }
        }
      }

    }

    /* If no boundary condition coefficients are available */

    else {

      for (cs_lnum_t i = 0; i < n_dist; i++) {
        cs_lnum_t f_id = dist_loc[i];
        cs_lnum_t c_id = b_face_cells[f_id];
        for (cs_lnum_t j = 0; j < dim; j++)
          distant_var[i*dim+j] = f->val[c_id*dim+j];
      }

    }

  }

  ple_locator_exchange_point_var(locator,
                                 distant_var,
                                 local_var,
                                 NULL,               /* faces indirection */
                                 sizeof(cs_real_t),
                                 f->dim,
                                 0);

  /* Now set boundary condition values */

  for (cs_lnum_t j = 0; j < dim; j++) {

    cs_real_t  *_rcodcl = rcodcl + (var_id+j)*n_b_faces;

    for (cs_lnum_t i = 0; i < n_faces; i++) {
      const cs_lnum_t f_id = (faces != NULL) ? faces[i] : i;
      _rcodcl[f_id] = local_var[i*dim + j];
    }

  }

  BFT_FREE(local_var);
  BFT_FREE(distant_var);

  /* Compute initial balance if applicable */

  if (normalize > 0) {

    _inlet_sum(var_id,
               f,
               cs_glob_mesh,
               cs_glob_mesh_quantities,
               normalize,
               n_faces,
               faces,
               balance_w,
               nvarcl,
               rcodcl,
               inlet_sum_1);

    for (cs_lnum_t j = 0; j < dim; j++) {

      const cs_real_t f_mult = (fabs(inlet_sum_1[j]) > 1.e-24) ?
                               inlet_sum_0[j] / inlet_sum_1[j] : 1.;

      cs_real_t  *_rcodcl = rcodcl + (var_id+j)*n_b_faces;

      for (cs_lnum_t i = 0; i < n_faces; i++) {
        const cs_lnum_t f_id = (faces != NULL) ? faces[i] : i;
        _rcodcl[f_id] *= f_mult;
      }

    }

  }
}

/*----------------------------------------------------------------------------
 * Create the boundary conditions type array bc_type
 *----------------------------------------------------------------------------*/

void
cs_boundary_conditions_type_create(void)
{
  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;

  BFT_MALLOC(_bc_type, n_b_faces, int);

  for (cs_lnum_t ii = 0 ; ii < n_b_faces ; ii++)
    _bc_type[ii] = 0;

  cs_glob_bc_type = _bc_type;
}

/*----------------------------------------------------------------------------
 * Free the boundary conditions type array _bc_type
 *----------------------------------------------------------------------------*/

void
cs_boundary_conditions_type_free(void)
{
  BFT_FREE(_bc_type);
}

/*----------------------------------------------------------------------------*/
/*! \brief Set convective oulet boundary condition for a scalar
 *
 * Parameters:
 * \param[out]    coefa         explicit BC coefficient for gradients
 * \param[out]    cofaf         explicit BC coefficient for diffusive flux
 * \param[out]    coefb         implicit BC coefficient for gradients
 * \param[out]    cofbf         implicit BC coefficient for diffusive flux
 * \param[in]     pimp          Flux value to impose
 * \param[in]     cfl           Local Courant number used to convect
 * \param[in]     hint          Internal exchange coefficient
 */
/*----------------------------------------------------------------------------*/

void cs_boundary_conditions_set_convective_outlet_scalar(cs_real_t *coefa ,
                                                         cs_real_t *cofaf,
                                                         cs_real_t *coefb,
                                                         cs_real_t *cofbf,
                                                         cs_real_t  pimp,
                                                         cs_real_t  cfl,
                                                         cs_real_t  hint)
{
  /* Gradient BCs */
  *coefb = cfl / (1.0 + cfl);
  *coefa = (1.0 - *coefb) * pimp;

  /* Flux BCs */
  *cofaf = - hint * *coefa;
  *cofbf =   hint * (1.0 - *coefb);

  return;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
