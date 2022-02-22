/*============================================================================
 * \file measures set management.
 *============================================================================*/

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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <float.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem_usage.h"
#include "bft_error.h"
#include "bft_printf.h"
#include "bft_mem.h"

#include "fvm_nodal.h"
#include "fvm_point_location.h"

#include "cs_base.h"
#include "cs_selector.h"
#include "cs_mesh_connect.h"
#include "cs_post.h"
#include "cs_prototypes.h"

#include "cs_log.h"
#include "cs_map.h"
#include "cs_parall.h"
#include "cs_mesh_location.h"

#include "cs_mesh.h"
#include "cs_measures_util.h"

#include "cs_field.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Measures set definitions */

static int  _n_measures_sets = 0;
static int  _n_measures_sets_max = 0;
static cs_measures_set_t  *_measures_sets = NULL;
static cs_map_name_to_id_t  *_measures_sets_map = NULL;

static int  _n_grids = 0;
static int  _n_grids_max = 0;
static cs_interpol_grid_t  *_grids = NULL;
static cs_map_name_to_id_t  *_grids_map = NULL;

/*============================================================================
 * Local Type Definitions
 *============================================================================*/

typedef struct
{
  double val;
  int    rank;
} _cs_base_mpi_double_int_t;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create mesh <-> interpolation grid connectivity
 *
 * parameters:
 *   ig             <-- pointer to the interpolation grid structure
 *----------------------------------------------------------------------------*/

static void
_mesh_interpol_create_connect(cs_interpol_grid_t   *ig)
{
  cs_lnum_t  ii;
  cs_lnum_t nb_points = ig->nb_points;
  cs_real_t *coords = ig->coords;

  cs_lnum_t *location = NULL;
  float *distance = NULL;
  fvm_nodal_t *nodal_mesh = NULL;
  const cs_mesh_t *mesh = cs_glob_mesh;

#if defined(HAVE_MPI)
  _cs_base_mpi_double_int_t val_min, val_in;
#endif

  nodal_mesh = cs_mesh_connect_cells_to_nodal(mesh,
                                              "temporary",
                                              false,
                                              mesh->n_cells,
                                              NULL);

  BFT_MALLOC(location, nb_points, cs_lnum_t);
  BFT_MALLOC(distance, nb_points, float);

#   pragma omp parallel for
  for (ii = 0; ii < nb_points; ii++) {
    location[ii] = -1;
    distance[ii] = -1.0;
  }

  fvm_point_location_nodal(nodal_mesh,
                           0,
                           0.1,
                           0,
                           nb_points,
                           NULL,
                           coords,
                           location,
                           distance);


#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    for (ii = 0; ii < nb_points; ii++) {
      if (location[ii] > 0)
        val_in.val = distance[ii];
      else
        val_in.val = DBL_MAX;

      val_in.rank = cs_glob_rank_id;

      MPI_Reduce(&val_in, &val_min, 1, MPI_DOUBLE_INT, MPI_MINLOC, 0,
                  cs_glob_mpi_comm);
      MPI_Bcast(&val_min.rank, 1, MPI_INT, 0, cs_glob_mpi_comm);
      MPI_Bcast(&location[ii], 1, MPI_INT, val_min.rank,
                cs_glob_mpi_comm);

      ig->rank_connect[ii] = val_min.rank;
    }
  }
#endif

  /* copy location and switch to 0-based */
#   pragma omp parallel for
  for (ii = 0; ii < nb_points; ii++) {
    ig->cell_connect[ii] = location[ii] - 1;
  }

  nodal_mesh = fvm_nodal_destroy(nodal_mesh);
  BFT_FREE(location);
  BFT_FREE(distance);
}

/*----------------------------------------------------------------------------
 * Interpolate mesh field on interpol grid structure.
 *
 * parameters:
 *   ig                   <-- pointer to the interpolation grid structure
 *   values_to_interpol   <-- field on mesh (size = n_cells)
 *   interpolated_values  --> interpolated values on the interpolation grid
 *                            structure (size = ig->nb_point)
 *----------------------------------------------------------------------------*/

void
cs_interpol_field_on_grid(cs_interpol_grid_t         *ig,
                          const cs_real_t            *values_to_interpol,
                          cs_real_t                  *interpoled_values)
{
  cs_lnum_t ii, jj;
  int ms_dim = 1;
  cs_lnum_t nb_points = ig->nb_points;
  cs_mesh_t *mesh = cs_glob_mesh;

#   pragma omp parallel for private(jj)
  for (ii = 0; ii < nb_points; ii++) {
    if (ig->cell_connect[ii] > -1 && ig->cell_connect[ii] < mesh->n_cells)
      for (jj = 0; jj < ms_dim; jj++)
        interpoled_values[ii*ms_dim +jj] =
          values_to_interpol[(ig->cell_connect[ii])*ms_dim + jj];
  }

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1)
    for (ii = 0; ii < nb_points; ii++)
      for (jj = 0; jj < ms_dim; jj++)
        MPI_Bcast(&(interpoled_values[ii*ms_dim + jj]), 1,
                  CS_MPI_REAL, ig->rank_connect[ii],
                  cs_glob_mpi_comm);

#endif
}

/*----------------------------------------------------------------------------
 * Compute a Cressman interpolation on the global mesh.
 *
 * parameters:
 *   ms                   <-- pointer to the measures set structure
 *                            (values to interpolate)
 *   interpolated_values  --> interpolated values on the global mesh
 *                            (size = n_cells or nb_faces)
 *   id_type              <-- parameter:
 *                              1: interpolation on volumes
 *                              2: interpolation on boundary faces
 *----------------------------------------------------------------------------*/

void
cs_cressman_interpol(cs_measures_set_t         *ms,
                     cs_real_t                 *interpolated_values,
                     int                        id_type)
{
  cs_lnum_t n_elts = 0;
  cs_real_t *xyz_cen = NULL;
  const cs_mesh_t *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t *mesh_quantities = cs_glob_mesh_quantities;

  if (id_type == 1) {
    n_elts = mesh->n_cells;
    xyz_cen = mesh_quantities->cell_cen;
  }
  else if (id_type == 2) {
    n_elts = mesh->n_b_faces;
    xyz_cen = mesh_quantities->b_face_cog;
  }
  assert(xyz_cen != NULL);

# pragma omp parallel for
  for (cs_lnum_t ii = 0; ii < n_elts; ii++) {
    cs_real_t total_weight = 0.;
    cs_real_t interpolated_value = 0.;
    for (cs_lnum_t jj = 0; jj < ms->nb_measures; jj++) {
      if (ms->is_cressman[jj] == 1) {
        cs_real_t dist_x = (xyz_cen[ii*3   ] - ms->coords[jj*3   ])
                          *ms->inf_radius[jj*3   ];
        cs_real_t dist_y = (xyz_cen[ii*3 +1] - ms->coords[jj*3 +1])
                          *ms->inf_radius[jj*3 +1];
        cs_real_t dist_z = (xyz_cen[ii*3 +2] - ms->coords[jj*3 +2])
                          *ms->inf_radius[jj*3 +2];

        cs_real_t r2 = dist_x*dist_x + dist_y*dist_y + dist_z*dist_z;

        cs_real_t weight = 0.;
        if (r2/4. <= 700.)
          weight = exp(-r2/4.);

        total_weight += weight;
        interpolated_value += (ms->measures[jj])*weight;
      }
    }

    if (total_weight > 0.)
      interpolated_values[ii] = interpolated_value/total_weight;
    else
      interpolated_values[ii] = 0.;
  }

}

/*----------------------------------------------------------------------------
 * Create an interpolation grid descriptor.
 *
 * For measures set with a dimension greater than 1, components are interleaved.
 *
 * parameters:
 *   name        <-- grid name
 *
 * returns:
 *   pointer to new interpolation grid.
 *
 *----------------------------------------------------------------------------*/

cs_interpol_grid_t *
cs_interpol_grid_create(const char   *name)
{
  bool reall = true;
  int grid_id = -1;

  const char *addr_0 = NULL, *addr_1 = NULL;

  cs_interpol_grid_t *ig =  NULL;

  /* Initialize if necessary */

  if (_grids_map == NULL)
    _grids_map = cs_map_name_to_id_create();
  else
    addr_0 = cs_map_name_to_id_reverse(_grids_map, 0);

  if (strlen(name) == 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Defining a interpolation grid requires a name."));

  /* Find or insert entry in map */

  grid_id = cs_map_name_to_id(_grids_map, name);

  /* Move name pointers of previous grid if necessary
     (i.e. reallocation of map names array) */

  addr_1 = cs_map_name_to_id_reverse(_grids_map, 0);

  if (addr_1 != addr_0) {
    int i;
    ptrdiff_t addr_shift = addr_1 - addr_0;
    for (i = 0; i < grid_id; i++)
      (_grids + i)->name += addr_shift;
  }

  if (grid_id == _n_grids) {
    _n_grids = grid_id + 1;
    reall = false;
  }

  /* Reallocate grids if necessary */

  if (_n_grids > _n_grids_max) {
    if (_n_grids_max == 0)
      _n_grids_max = 8;
    else
      _n_grids_max *= 2;
    BFT_REALLOC(_grids, _n_grids_max, cs_interpol_grid_t);
  }

  /* Assign grid */

  ig = _grids + grid_id;

  ig->name = cs_map_name_to_id_reverse(_grids_map, grid_id);

  ig->id = grid_id;

  ig->nb_points = 0;

  if (!reall) {
    ig->coords = NULL;
    ig->cell_connect = NULL;
    ig->rank_connect = NULL;
  }

  else {
    BFT_FREE(ig->coords);
    if (ig->is_connect) {
      BFT_FREE(ig->cell_connect);
#if defined(HAVE_MPI)
      if (cs_glob_n_ranks > 1)
        BFT_FREE(ig->rank_connect);
#endif
    }
  }
  ig->is_connect = false;

  return ig;
}

/*----------------------------------------------------------------------------
 * Create interpolation grid.
 *
 * parameters:
 *   ig        <-- pointer to an interpol grid structure
 *   nb_points <-- number of considered points
 *   coord     <-- coordonates of considered points
 *----------------------------------------------------------------------------*/

void
cs_interpol_grid_init(cs_interpol_grid_t    *ig,
                      const cs_lnum_t        nb_points,
                      const cs_real_t       *coords)
{
  cs_lnum_t ii;
  BFT_MALLOC(ig->cell_connect, nb_points, cs_lnum_t);
#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1)
    BFT_MALLOC(ig->rank_connect, nb_points, int);
#endif
  BFT_MALLOC(ig->coords, 3*nb_points, cs_real_t);

#   pragma omp parallel for
  for (ii = 0; ii < 3*nb_points; ii++) {
    ig->coords[ii] = coords[ii];
  }

  ig->nb_points = nb_points;

  _mesh_interpol_create_connect(ig);

  ig->is_connect = true;
}

/*----------------------------------------------------------------------------
 * Create a measures set descriptor.
 *
 * For measures set with a dimension greater than 1, components are interleaved.
 *
 * parameters:
 *   name        <-- measures set name
 *   type_flag   <-- mask of field property and category values (not used yet)
 *   dim         <-- measure set dimension (number of components)
 *   interleaved <-- if dim > 1, indicate if field is interleaved
 *
 * returns:
 *   pointer to new measures set.
 *
 *----------------------------------------------------------------------------*/

cs_measures_set_t *
cs_measures_set_create(const char   *name,
                       int           type_flag,
                       int           dim,
                       bool          interleaved)
{
  bool reall = true;
  int measures_set_id = -1;

  const char *addr_0 = NULL, *addr_1 = NULL;

  cs_measures_set_t *ms =  NULL;

  /* Initialize if necessary */

  if (_measures_sets_map == NULL)
    _measures_sets_map = cs_map_name_to_id_create();
  else
    addr_0 = cs_map_name_to_id_reverse(_measures_sets_map, 0);

  if (strlen(name) == 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Defining a measure set requires a name."));

  /* Find or insert entry in map */

  measures_set_id = cs_map_name_to_id(_measures_sets_map, name);

  /* Move name pointers of previous measure set if necessary
     (i.e. reallocation of map names array) */

  addr_1 = cs_map_name_to_id_reverse(_measures_sets_map, 0);

  if (addr_1 != addr_0) {
    int i;
    ptrdiff_t addr_shift = addr_1 - addr_0;
    for (i = 0; i < measures_set_id; i++)
      (_measures_sets + i)->name += addr_shift;
  }

  if (measures_set_id == _n_measures_sets) {
    _n_measures_sets = measures_set_id + 1;
    reall = false;
  }

  /* Reallocate measures sets if necessary */

  if (_n_measures_sets > _n_measures_sets_max) {
    if (_n_measures_sets_max == 0)
      _n_measures_sets_max = 8;
    else
      _n_measures_sets_max *= 2;
    BFT_REALLOC(_measures_sets, _n_measures_sets_max, cs_measures_set_t);
  }

  /* Assign measure set */

  ms = _measures_sets + measures_set_id;

  ms->name = cs_map_name_to_id_reverse(_measures_sets_map, measures_set_id);

  ms->id = measures_set_id;
  ms->type = type_flag;
  ms->dim = dim;

  if (ms->dim > 1)
    ms->interleaved = interleaved;
  else
    ms->interleaved = true;

  ms->nb_measures = 0;
  ms->nb_measures_max = 0;

  if (!reall) {
    ms->coords = NULL;
    ms->measures = NULL;
    ms->is_cressman = NULL;
    ms->is_interpol = NULL;
    ms->inf_radius = NULL;
    ms->comp_ids = NULL;
  }
  else {
    BFT_FREE(ms->coords);
    BFT_FREE(ms->measures);
    BFT_FREE(ms->is_cressman);
    BFT_FREE(ms->is_interpol);
    BFT_FREE(ms->inf_radius);
    BFT_FREE(ms->comp_ids);
  }

  return ms;
}

/*----------------------------------------------------------------------------
 * (re)Allocate and fill in a measures set structure with an array of measures.
 *
 * parameters:
 *   ms               <-- pointer to the measures set
 *   nb_measures      <-- number of measures
 *   is_cressman      <-- for each measure cressman interpolation is:
 *                          0: not used
 *                          1: used
 *   is_interpol      <-- for the interpolation on mesh, each measure is
 *                          0: not taken into account
 *                          1: taken into account
 *   measures_coords  <-- measures spaces coordonates
 *   measures         <-- measures values (associated to coordinates)
 *   influence_radius <-- influence radius for interpolation (xyz interleaved)
 *----------------------------------------------------------------------------*/

void
cs_measures_set_map_values(cs_measures_set_t       *ms,
                           const cs_lnum_t          nb_measures,
                           const int               *is_cressman,
                           const int               *is_interpol,
                           const cs_real_t         *measures_coords,
                           const cs_real_t         *measures,
                           const cs_real_t         *influence_radius)
{
  cs_lnum_t  ii;
  int dim = ms->dim;

  if (nb_measures != ms->nb_measures) {
    BFT_REALLOC(ms->measures, nb_measures*dim, cs_real_t);
    BFT_REALLOC(ms->inf_radius, nb_measures*3, cs_real_t);
    BFT_REALLOC(ms->coords, nb_measures*3, cs_real_t);
    BFT_REALLOC(ms->is_cressman, nb_measures, int);
    BFT_REALLOC(ms->is_interpol, nb_measures, int);
    ms->nb_measures = nb_measures;
  }

  if (dim == 1) {
#   pragma omp parallel for
    for (ii = 0; ii < nb_measures; ii++)
      ms->measures[ii] = measures[ii];
  }
  else {
    if (ms->interleaved) {
      cs_lnum_t jj;
#   pragma omp parallel for private(jj)
      for (ii = 0; ii < nb_measures; ii++) {
        for (jj = 0; jj < dim; jj++)
          ms->measures[ii*dim + jj] = measures[ii*dim + jj];
      }
    }
    else {
      cs_lnum_t jj;
#   pragma omp parallel for private(jj)
      for (ii = 0; ii < nb_measures; ii++) {
        for (jj = 0; jj < dim; jj++)
          ms->measures[ii*dim + jj] = measures[jj*nb_measures + ii];
      }
    }
  }
#   pragma omp parallel for
  for (ii = 0; ii < nb_measures; ii++) {
    ms->is_interpol[ii] = is_interpol[ii];
    ms->is_cressman[ii] = is_cressman[ii];
  }

  {
    cs_lnum_t jj;
#   pragma omp parallel for private(jj)
    for (ii = 0; ii < nb_measures; ii++)
      for (jj = 0; jj < 3; jj++) {
        ms->coords[ii*3 + jj] = measures_coords[ii*3 + jj];
        ms->inf_radius[ii*3 +jj] = influence_radius[ii*3 + jj];
      }
  }
}

/*----------------------------------------------------------------------------
 * Add new measures to an existing measures set (already declared and
 * allocated).
 *
 * parameters:
 *   ms               <-- pointer to the existing measures set
 *   nb_measures      <-- number of new measures
 *   is_cressman      <-- for each new measure cressman interpolation is:
 *                          0: not used
 *                          1: used
 *   is_interpol      <-- for the interpolation on mesh, each new measure is
 *                          0: not taken into account
 *                          1: taken into account
 *   measures_coords  <-- new measures spaces coordonates
 *   measures         <-- new measures values (associated to coordonates)
 *   influence_radius <-- influence radius for interpolation (xyz interleaved)
 *----------------------------------------------------------------------------*/

void
cs_measures_set_add_values(cs_measures_set_t       *ms,
                           const cs_lnum_t          nb_measures,
                           const int               *is_cressman,
                           const int               *is_interpol,
                           const cs_real_t         *measures_coords,
                           const cs_real_t         *measures,
                           const cs_real_t         *influence_radius)
{
  cs_lnum_t  ii;
  int dim = ms->dim;
  if (ms->nb_measures + nb_measures > ms->nb_measures_max) {
    ms->nb_measures_max = 2*(ms->nb_measures + nb_measures);

    BFT_REALLOC(ms->measures, ms->nb_measures_max*dim, cs_real_t);
    BFT_REALLOC(ms->coords, ms->nb_measures_max*3, cs_real_t);
    BFT_REALLOC(ms->is_cressman, ms->nb_measures_max, int);
    BFT_REALLOC(ms->is_interpol, ms->nb_measures_max, int);
  }

  if (dim == 1) {
#   pragma omp parallel for
    for (ii = 0; ii < nb_measures; ii++)
      ms->measures[ms->nb_measures + ii] = measures[ii];
  }
  else {
    if (ms->interleaved) {
      cs_lnum_t jj;
#   pragma omp parallel for private(jj)
      for (ii = 0; ii < nb_measures; ii++) {
        for (jj = 0; jj < dim; jj++)
          ms->measures[(ii + ms->nb_measures)*dim + jj] = measures[ii*dim + jj];
      }
    }
    else {
      cs_lnum_t jj;
#   pragma omp parallel for private(jj)
      for (ii = 0; ii < nb_measures; ii++) {
        for (jj = 0; jj < dim; jj++)
          ms->measures[ii*dim + jj] = measures[ii*nb_measures + jj];
      }
    }
  }
#   pragma omp parallel for
  for (ii = 0; ii < nb_measures; ii++) {
    ms->is_interpol[ii + ms->nb_measures] = is_interpol[ii];
    ms->is_cressman[ii + ms->nb_measures] = is_cressman[ii];
  }

  {
    cs_lnum_t jj;
#   pragma omp parallel for private(jj)
    for (ii = 0; ii < nb_measures; ii++)
      for (jj = 0; jj < 3; jj++) {
        ms->coords[(ii + ms->nb_measures)*3 + jj] = measures_coords[ii*3 + jj];
        ms->inf_radius[(ii + ms->nb_measures)*3 + jj] = influence_radius[ii*3 + jj];
      }
  }

  ms->nb_measures += nb_measures;
}

/*----------------------------------------------------------------------------
 * Return a pointer to a measures set based on its id.
 *
 * This function requires that a measures set of the given id is defined.
 *
 * parameters:
 *  id <-- measures set id
 *
 * return:
 *   pointer to the measures set structure
 *
 *----------------------------------------------------------------------------*/

cs_measures_set_t  *
cs_measures_set_by_id(int  id)
{
  if (id > -1 && id < _n_measures_sets)
    return _measures_sets + id;
  else {
    bft_error(__FILE__, __LINE__, 0,
              _("Measure set with id %d is not defined."), id);
    return NULL;
  }
}

/*----------------------------------------------------------------------------
 * Return a pointer to a grid based on its id.
 *
 * This function requires that a grid of the given id is defined.
 *
 * parameters:
 *  id <-- grid id
 *
 * return:
 *   pointer to the grid structure
 *
 *----------------------------------------------------------------------------*/

cs_interpol_grid_t  *
cs_interpol_grid_by_id(int  id)
{
  if (id > -1 && id < _n_grids)
    return _grids + id;
  else {
    bft_error(__FILE__, __LINE__, 0,
              _("Interpol grid with id %d is not defined."), id);
    return NULL;
  }
}

/*----------------------------------------------------------------------------
 * Return a pointer to a measure set based on its name.
 *
 * This function requires that a measure set of the given name is defined.
 *
 * parameters:
 * name <-- measure set name
 *
 * return:
 *   pointer to the measures set structure
 *----------------------------------------------------------------------------*/

cs_measures_set_t  *
cs_measures_set_by_name(const char  *name)
{
  int id = cs_map_name_to_id_try(_measures_sets_map, name);

  if (id > -1)
    return _measures_sets + id;
  else {
    bft_error(__FILE__, __LINE__, 0,
              _("Measure set \"%s\" is not defined."), name);
    return NULL;
  }
}

/*----------------------------------------------------------------------------
 * Return a pointer to a grid based on its name.
 *
 * This function requires that a grid of the given name is defined.
 *
 * parameters:
 * name <-- grid name
 *
 * return:
 *   pointer to the grid structure
 *----------------------------------------------------------------------------*/

cs_interpol_grid_t  *
cs_interpol_grid_by_name(const char  *name)
{
  int id = cs_map_name_to_id_try(_grids_map, name);

  if (id > -1)
    return _grids + id;
  else {
    bft_error(__FILE__, __LINE__, 0,
              _("Interpol grid \"%s\" is not defined."), name);
    return NULL;
  }
}


/*----------------------------------------------------------------------------
 * Destroy all defined measures sets.
 *----------------------------------------------------------------------------*/

void
cs_measures_sets_destroy(void)
{
  int i;

  for (i = 0; i < _n_measures_sets; i++) {
    cs_measures_set_t  *ms = _measures_sets + i;
    BFT_FREE(ms->measures);
    BFT_FREE(ms->coords);
    BFT_FREE(ms->is_interpol);
    BFT_FREE(ms->is_cressman);
    BFT_FREE(ms->inf_radius);
    BFT_FREE(ms->comp_ids);
  }

  BFT_FREE(_measures_sets);

  cs_map_name_to_id_destroy(&_measures_sets_map);

  _n_measures_sets = 0;
  _n_measures_sets_max = 0;
}

/*----------------------------------------------------------------------------
 * Destroy all defined grids.
 *----------------------------------------------------------------------------*/

void
cs_interpol_grids_destroy(void)
{
  int i;

  for (i = 0; i < _n_grids; i++) {
    cs_interpol_grid_t  *ig = _grids + i;
    BFT_FREE(ig->coords);
    BFT_FREE(ig->cell_connect);
    if (cs_glob_n_ranks > 1)
      BFT_FREE(ig->rank_connect);
  }

  BFT_FREE(_grids);

  cs_map_name_to_id_destroy(&_grids_map);

  _n_grids = 0;
  _n_grids_max = 0;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define a measures set.
 *
 * Fortran interface; use mestcr;
 *
 * subroutine mestcr (name, idim, ilved, imeset)
 * *****************
 *
 * character*       name        : <-- : Measure set name
 * integer          idim        : <-- : Measures set dimension
 * integer          ilved       : <-- : 0: not intereaved; 1: interleaved
 * integer          imesset     : --> : id of defined measures set
 *----------------------------------------------------------------------------*/

void CS_PROCF(mestcr, MESTCR)
(
 const char   *name,
 const int    *idim,
 const int    *ilved,
 int          *imeset
)
{
  int type_flag = 0;
  bool interleaved = (*ilved == 0) ? false : true;
  cs_measures_set_t *ms = NULL;

  ms = cs_measures_set_create(name,
                              type_flag,
                              *idim,
                              interleaved);

  *imeset = ms->id;
}

/*----------------------------------------------------------------------------
 * Define a grid.
 *
 * Fortran interface
 *
 * subroutine gridcr (name, igrid)
 * *****************
 *
 * character*       name        : <-- : Measure set name
 * integer          igrid       : --> : id of defined grid
 *----------------------------------------------------------------------------*/

void CS_PROCF(gridcr, GRIDCR)
(
 const char     *name,
 int            *igrid
)
{
  cs_interpol_grid_t *ig = NULL;

  ig = cs_interpol_grid_create(name);

  *igrid = ig->id;
}

/*----------------------------------------------------------------------------
 * (re)Allocate and map values to a measure set.
 *
 * Fortran interface
 *
 * subroutine mesmap (imeset, inbmes, meset, coords, cressm, interp)
 * *****************
 *
 * integer          imeset      : <-- : Measures set id
 * integer          inbmes      : <-- : Number of measures
 * cs_real_t*       meset       : <-- : Pointer to measures values array
 * cs_real_t*       coords      : <-- : Pointer to measures coordonates array
 * integer*         cressm      : <-- : Pointer to Cressman interpolation flag
 * integer*         interp      : <-- : Pointer to interpolation flag
 * integer*         infrad      : <-- : Influence radius for interpolation
 *----------------------------------------------------------------------------*/

void CS_PROCF(mesmap, MESMAP)
(
 const int         *imeset,
 const int         *inbmes,
 const cs_real_t   *meset,
 const cs_real_t   *coords,
 const int         *cressm,
 const int         *interp,
 const cs_real_t   *infrad
)
{

  cs_measures_set_t *ms = cs_measures_set_by_id(*imeset);

  cs_measures_set_map_values(ms,
                             *inbmes,
                             cressm,
                             interp,
                             coords,
                             meset,
                             infrad);

}

/*----------------------------------------------------------------------------
 * Map a grid grid.
 *
 * Fortran interface
 *
 * subroutine gridmap (name, lname, igrid)
 * *****************
 *
 * integer          igrid       : <-- : Measures set id
 * integer          inpts       : <-- : Number of measures
 * cs_real_t*       coords      : <-- : Pointer to measures coordonates array
 *----------------------------------------------------------------------------*/

void CS_PROCF(grimap, GRIMAP)
(
 const int         *igrid,
 const int         *inpts,
 const cs_real_t   *coords
)
{
  cs_interpol_grid_t *ig = cs_interpol_grid_by_id(*igrid);

  cs_interpol_grid_init(ig,
                        *inpts,
                        coords);
}

/*----------------------------------------------------------------------------
 * Add values to a measure set.
 *
 * Fortran interface
 *
 * subroutine mesadd (imeset, inbmes, meset, coords, cressm, interp)
 * *****************
 *
 * integer          imeset      : <-- : Measures set id
 * integer          inbmes      : <-- : Number of measures to add
 * cs_real_t*       meset       : <-- : Pointer to measures values array
 * cs_real_t*       coords      : <-- : Pointer to measures coordonates array
 * integer*         cressm      : <-- : Pointer to Cressman interpolation flag
 * integer*         interp      : <-- : Pointer to interpolation flag
 * integer*         infrad      : <-- : Influence radius for interpolation
 *----------------------------------------------------------------------------*/

void CS_PROCF(mesadd, MESADD)
(
 const int         *imeset,
 const int         *inbmes,
 const cs_real_t   *meset,
 const cs_real_t   *coords,
 const int         *cressm,
 const int         *interp,
 const cs_real_t   *infrad
)
{
  cs_measures_set_t *ms = cs_measures_set_by_id(*imeset);

  cs_measures_set_add_values(ms,
                             *inbmes,
                             cressm,
                             interp,
                             coords,
                             meset,
                             infrad);

}

/*----------------------------------------------------------------------------
 * Compute a Cressman interpolation on the global mesh.
 *
 * Fortran interface
 *
 * subroutine mscrss (imeset, type, pldval)
 * *****************
 *
 * integer          imeset      : <-- : Measures set id
 * integer          type        : <-- : Parameter:
 *                                        1: interpolation on volumes
 *                                        2: interpolation on boundary faces
 * cs_real_t*       pldval      : --> : Interpolated values on the global mesh
 *----------------------------------------------------------------------------*/

void CS_PROCF(mscrss, MSCRSS)
(
 const int         *imeset,
 const int         *type,
 cs_real_t         *pldval
)
{
  cs_measures_set_t *ms = cs_measures_set_by_id(*imeset);

  cs_cressman_interpol(ms,
                       pldval,
                       *type);
}

/*----------------------------------------------------------------------------
 * Interpolate calculed field on a grid.
 *
 * Fortran interface
 *
 * subroutine gripol (igrid, inval, pldval)
 * *****************
 *
 * integer          igrid       : <-- : Measures set id
 * cs_real_t*       inval       : <-- : Values to interpolate
 * cs_real_t*       pldval      : --> : Interpolated values on the grid
 *----------------------------------------------------------------------------*/

void CS_PROCF(gripol, GRIPOL)
(
 const int         *igrid,
 const cs_real_t   *inval,
 cs_real_t         *pldval
)
{
  cs_interpol_grid_t *ig = cs_interpol_grid_by_id(*igrid);

  cs_interpol_field_on_grid(ig,
                            inval,
                            pldval);
}

/*----------------------------------------------------------------------------
 * Destroy measures sets.
 *
 * Fortran interface
 *
 * subroutine mestde (void)
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF(mestde, MESTDE)(void)
{
  cs_measures_sets_destroy();
}

/*----------------------------------------------------------------------------
 * Destroy grids.
 *
 * Fortran interface
 *
 * subroutine grides (void)
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF(grides, GRIDES)(void)
{
  cs_interpol_grids_destroy();
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
