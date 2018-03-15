/*============================================================================
 * Boundary condition handling.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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
#include "cs_gui_util.h"
#include "cs_field.h"
#include "cs_field_operator.h"
#include "cs_field_pointer.h"
#include "cs_flag_check.h"
#include "cs_halo.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_connect.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_physical_model.h"
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

/*----------------------------------------------------------------------------
 * Local Structure Definitions
 *----------------------------------------------------------------------------*/

/*============================================================================
 *  Global variables
 *============================================================================*/

static int *_bc_type;

const int *cs_glob_bc_type = NULL;

static int *_bc_face_zone;

const int *cs_glob_bc_face_zone = NULL;

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
                                    int                        nvar,
                                    cs_real_t                 *rcodcl);

void
cs_f_boundary_conditions_get_pointers(int  **itypfb,
                                      int  **izfppp);

/*============================================================================
 * Private function definitions
 *============================================================================*/

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
 *   nvar            <-- number of variables requiring BC's
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
           int                          nvar,
           cs_real_t                    rcodcl[],
           cs_real_t                    inlet_sum[])
{
  CS_UNUSED(nvar);

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
 *   nvar            <-- number of variables requiring BC's
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
                                    int                        nvar,
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
                                    nvar,
                                    rcodcl);

  BFT_FREE(_faces);
}

/*----------------------------------------------------------------------------
 * Get pointer to cs_glob_bc_type and cs_glob_bc_face_zone
 *----------------------------------------------------------------------------*/

void
cs_f_boundary_conditions_get_pointers(int **itypfb,
                                      int **izfppp)
{
  *itypfb = _bc_type;
  *izfppp = _bc_face_zone;
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
 * \param[in]  bc_flag     array of BC type ids
 * \param[in]  type_name   name of attribute in error, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_error(const int       *bc_flag,
                             const char      *type_name)
{
  const char type_name_default[] = N_("boundary condition type");
  const char *_type_name = type_name_default;

  if (type_name != NULL)
    _type_name = type_name;

  /* Check for faces with problems;
     bc_flag[] used to determine if we have an error */

  int have_errors
    = cs_flag_check(_("face with boundary condition definition error"),
                    _type_name,
                    _("BC type"),
                    _("Faces with B.C. error"),
                    _("Faces with valid B.C.'s"),
                    CS_MESH_LOCATION_BOUNDARY_FACES,
                    1, /* min_flag */
                    bc_flag);

  if (have_errors)
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
 * \param[in]       nvar             number of variables requiring BC's
 * \param[in, out]  rcodcl           boundary condition values
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_mapped_set(const cs_field_t          *f,
                                  ple_locator_t             *locator,
                                  cs_mesh_location_type_t    location_type,
                                  int                        normalize,
                                  int                        interpolate,
                                  cs_lnum_t                  n_faces,
                                  const cs_lnum_t           *faces,
                                  cs_real_t                 *balance_w,
                                  int                        nvar,
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
               nvar,
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

  if (location_type == CS_MESH_LOCATION_CELLS || interpolate) {
    /* FIXME: we cheat here with the constedness of the field relative to
       for a possible ghost values update, but having a finer control
       of when syncing is required would be preferable */
    cs_field_t *_f = cs_field_by_id(f->id);
    cs_field_interpolate(_f,
                         interpolation_type,
                         n_dist,
                         dist_loc,
                         (const cs_real_3_t *)dist_coords,
                         distant_var);
  }

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
               nvar,
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
 * Create the boundary conditions face type and face zone arrays
 *----------------------------------------------------------------------------*/

void
cs_boundary_conditions_create(void)
{
  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;

  /* boundary conditions type by boundary face */

  BFT_MALLOC(_bc_type, n_b_faces, int);
  for (cs_lnum_t ii = 0 ; ii < n_b_faces ; ii++) {
    _bc_type[ii] = 0;
  }
  cs_glob_bc_type = _bc_type;

  /* boundary conditions zone by boundary face */
  /* only for specific physics */

  if (   cs_glob_physical_model_flag[CS_PHYSICAL_MODEL_FLAG] > 0
      || cs_gui_file_is_loaded()) {
    BFT_MALLOC(_bc_face_zone, n_b_faces, int);
    for (cs_lnum_t ii = 0 ; ii < n_b_faces ; ii++) {
      _bc_face_zone[ii] = 0;
    }
    cs_glob_bc_face_zone = _bc_face_zone;
  }
}

/*----------------------------------------------------------------------------
 * Free the boundary conditions face type and face zone arrays
 *----------------------------------------------------------------------------*/

void
cs_boundary_conditions_free(void)
{
  BFT_FREE(_bc_type);

  if (   cs_glob_physical_model_flag[CS_PHYSICAL_MODEL_FLAG] > 0
      || cs_gui_file_is_loaded()) {
    BFT_FREE(_bc_face_zone);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set convective oulet boundary condition for a scalar.
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

void
cs_boundary_conditions_set_convective_outlet_scalar(cs_real_t *coefa ,
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
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set free oulet boundary condition for the pressure.
 *
 * \param[in]     visc_term          \f$ \divv (\mu \gradv \vect{u}) \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_set_free_outlet(cs_real_3_t *visc_term)
{
  const cs_real_3_t  *cell_cen
        = (const cs_real_3_t *)cs_glob_mesh_quantities->cell_cen;
  const cs_real_3_t *surfac =
    (const cs_real_3_t *)cs_glob_mesh_quantities->i_face_normal;
  const cs_real_3_t *surfbo =
    (const cs_real_3_t *)cs_glob_mesh_quantities->b_face_normal;
  const cs_real_t *surfbn = cs_glob_mesh_quantities->b_face_surf;
  const cs_real_t *b_dist = cs_glob_mesh_quantities->b_dist;
  const cs_real_3_t *b_face_cog = (const cs_real_3_t *)cs_glob_mesh_quantities->b_face_cog;
  const cs_lnum_2_t *i_face_cells =
    (const cs_lnum_2_t *)(cs_glob_mesh->i_face_cells);
  const cs_lnum_t *b_face_cells = cs_glob_mesh->b_face_cells;
  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;
  const cs_lnum_t n_i_faces = cs_glob_mesh->n_i_faces;
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_real_t *weight = cs_glob_mesh_quantities->weight;
  const int key_cal_opt_id = cs_field_key_id("var_cal_opt");
  cs_var_cal_opt_t vcoptp;
  cs_real_t *pimp, *hint, pfac, pbik;
  cs_lnum_t ii, jj, if1, if2, iel, isou;
  cs_lnum_t *face_ids_outlet;

  cs_field_t *f = CS_F_(p);

  BFT_MALLOC(face_ids_outlet, n_cells, cs_lnum_t);
  BFT_MALLOC(pimp, n_b_faces, cs_real_t);
  BFT_MALLOC(hint, n_b_faces, cs_real_t);

# pragma omp parallel for
  for (iel = 0 ; iel < n_cells ; iel++)
    face_ids_outlet[iel] = 0;

# pragma omp parallel for
  for (cs_lnum_t ifac = 0 ; ifac < n_b_faces ; ifac++) {
    iel = b_face_cells[ifac];
    pimp[ifac] = 0.;
    if (cs_glob_bc_type[ifac] == CS_FREE_OUTLET)
      face_ids_outlet[iel] = ifac;
  }

  /* Get the calculation option from the pressure field */
  cs_field_get_key_struct(CS_F_(p), key_cal_opt_id, &vcoptp);

  cs_real_6_t *dttens, *tensor_permea;
  cs_real_33_t visci;
  cs_real_t *permea, viscis, fikis, distfi, distbf;

  if (vcoptp.idften == 3 || vcoptp.idften == 6)
    dttens = (cs_real_6_t *)cs_field_by_name("dttens")->val;

  if (cs_glob_physical_model_flag[CS_GROUNDWATER] == 1)
    permea = cs_field_by_name("permeability")->val;

  /* Computation of hint */
# pragma omp parallel for
  for (cs_lnum_t ifac = 0 ; ifac < n_b_faces ; ifac++) {

    iel = b_face_cells[ifac];
    distbf = b_dist[ifac];

    /* if a flux dt.grad p (w/m2) is set in cs_user_boundary */
    if (vcoptp.idften == 1) {
      hint[ifac] = CS_F_(dt)->val[iel]/distbf;
      if (cs_glob_physical_model_flag[CS_GROUNDWATER] == 1)
        hint[ifac] = permea[iel]/distbf;
    } else if (vcoptp.idften == 3) {
      hint[ifac] = (dttens[iel][0]*cs_math_sq(surfbo[ifac][0])
                 +  dttens[iel][1]*cs_math_sq(surfbo[ifac][1])
                 +  dttens[iel][2]*cs_math_sq(surfbo[ifac][2]))
                 / (cs_math_sq(surfbn[ifac]) * distbf);
    } else if (vcoptp.idften == 6) {
      if (cs_glob_physical_model_flag[CS_GROUNDWATER] == -1) {
        visci[0][0] = dttens[iel][0];
        visci[1][1] = dttens[iel][1];
        visci[2][2] = dttens[iel][2];
        visci[1][0] = dttens[iel][3];
        visci[0][1] = dttens[iel][3];
        visci[1][2] = dttens[iel][4];
        visci[2][1] = dttens[iel][4];
        visci[0][2] = dttens[iel][5];
        visci[2][0] = dttens[iel][5];
      } else {
        visci[0][0] = tensor_permea[iel][0];
        visci[1][1] = tensor_permea[iel][1];
        visci[2][2] = tensor_permea[iel][2];
        visci[1][0] = tensor_permea[iel][3];
        visci[0][1] = tensor_permea[iel][3];
        visci[1][2] = tensor_permea[iel][4];
        visci[2][1] = tensor_permea[iel][4];
        visci[0][2] = tensor_permea[iel][5];
        visci[2][0] = tensor_permea[iel][5];
      }

      /* ||ki.s||^2 */
      viscis = cs_math_sq(visci[0][0]*surfbo[ifac][0]
                         +visci[1][0]*surfbo[ifac][1]
                         +visci[2][0]*surfbo[ifac][2])
             + cs_math_sq(visci[0][1]*surfbo[ifac][0]
                         +visci[1][1]*surfbo[ifac][1]
                         +visci[2][1]*surfbo[ifac][2])
             + cs_math_sq(visci[0][2]*surfbo[ifac][0]
                         +visci[1][2]*surfbo[ifac][1]
                         +visci[2][2]*surfbo[ifac][2]);

      fikis = ((b_face_cog[ifac][0]-cell_cen[iel][0])*visci[0][0]
              +(b_face_cog[ifac][1]-cell_cen[iel][1])*visci[0][1]
              +(b_face_cog[ifac][2]-cell_cen[iel][2])*visci[0][2])*surfbo[ifac][0]
            + ((b_face_cog[ifac][0]-cell_cen[iel][0])*visci[1][0]
              +(b_face_cog[ifac][1]-cell_cen[iel][1])*visci[1][1]
              +(b_face_cog[ifac][2]-cell_cen[iel][2])*visci[1][2])*surfbo[ifac][1]
            + ((b_face_cog[ifac][0]-cell_cen[iel][0])*visci[2][0]
              +(b_face_cog[ifac][1]-cell_cen[iel][1])*visci[2][1]
              +(b_face_cog[ifac][2]-cell_cen[iel][2])*visci[2][2])*surfbo[ifac][2];

      distfi = b_dist[ifac];

      /* take i" so that i"f= eps*||fi||*ki.n when j" is in cell rji */
      /* nb: eps = 0.1 must be consistent with vitens.f90 */
      fikis = CS_MAX(fikis, 0.1*cs_math_sq(viscis)*distfi);
      hint[ifac] = viscis/surfbn[ifac]/fikis;
    }

    if (cs_field_by_name_try("void_fraction") != NULL)
      hint[ifac] /= CS_F_(rho)->val[iel];
  }

  /* Free outlet : interior faces contribution */
# pragma omp parallel for private(isou)
  for (cs_lnum_t ifac = 0 ; ifac < n_i_faces ; ifac++) {
    ii = i_face_cells[ifac][0];
    jj = i_face_cells[ifac][1];
    if1 = face_ids_outlet[ii];
    if2 = face_ids_outlet[jj];
    pfac = weight[ifac]*f->val_pre[ii] + (1.-weight[ifac])*f->val_pre[jj];
    for (isou = 0 ; isou < 3 ; isou++) {
      if (if1 > 0) pimp[if1] -= pfac*surfac[ifac][isou]*surfbo[if1][isou]/cs_math_sq(surfbn[if1]);
      if (if2 > 0) pimp[if2] += pfac*surfac[ifac][isou]*surfbo[if2][isou]/cs_math_sq(surfbn[if2]);
    }
  }

  /* Free outlet : viscous terms contribution */
# pragma omp parallel for private(isou)
  for (iel = 0 ; iel < n_cells ; iel++) {
    if1 = face_ids_outlet[iel];
    if (if1 > 0) {
      for (isou = 0 ; isou < 3 ; isou++)
        pimp[if1] += visc_term[iel][isou]*surfbo[if1][isou]/cs_math_sq(surfbn[if1]);
    }
  }

  /* Free outlet : contribution of the other boundary faces */
  /* TODO : reconstruct f->val_pre */
# pragma omp parallel for private(isou)
  for (cs_lnum_t ifac = 0 ; ifac < n_b_faces ; ifac++) {
    if (cs_glob_bc_type[ifac] != CS_FREE_OUTLET) {
      iel = b_face_cells[ifac];
      if1 = face_ids_outlet[iel];
      pbik = f->bc_coeffs->a[ifac] + f->bc_coeffs->b[ifac]*f->val_pre[iel];
      if (if1 > 0) {
        for (isou = 0 ; isou < 3 ; isou++)
           pimp[if1] -= pbik*surfbo[ifac][isou]*surfbo[if1][isou]/cs_math_sq(surfbn[if1]);
      }
    }
  }

  /* Free outlet : boundary coefficient computation */
# pragma omp parallel for
  for (cs_lnum_t ifac = 0 ; ifac < n_b_faces ; ifac++) {
    if (cs_glob_bc_type[ifac] == CS_FREE_OUTLET) {
      /* Gradient BCs */
      f->bc_coeffs->a[ifac] = pimp[ifac];
      f->bc_coeffs->b[ifac] = 0.;

      /* Flux BCs */
      f->bc_coeffs->af[ifac] = - hint[ifac] * pimp[ifac];
      f->bc_coeffs->bf[ifac] =   hint[ifac];
    }
  }

  BFT_FREE(face_ids_outlet);
  BFT_FREE(pimp);
  BFT_FREE(hint);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
