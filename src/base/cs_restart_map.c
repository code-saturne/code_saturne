/*============================================================================
 * Checkpoint / restart extension to mapped meshes
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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "ple_locator.h"

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"
#include "cs_array.h"
#include "cs_assert.h"
#include "cs_base.h"
#include "cs_coupling.h"
#include "cs_io.h"
#include "cs_mesh.h"
#include "cs_mesh_connect.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_preprocessor_data.h"
#include "cs_timer.h"
#include "cs_timer_stats.h"
#include "fvm_nodal.h"
#include "fvm_interpolate.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_restart.h"
#include "cs_restart_map.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Local type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Apply mesh deformation to previous mesh for mapping if present ? */

static bool _apply_mesh_deformation = false;

/* Variables to be saved between definition of mapping and
   build of that mapping */

static char *_mesh_input_path = NULL;
static float _tolerance[2] = {0, 0.1};

/* Save previous section read function */
static cs_restart_read_section_t   *_read_section_f  = NULL;

/* PLE locators for main mesh locations (not all used) */
bool _need_locator[4] = {true, false, false, false};
static  ple_locator_t  *_locator[4] = {NULL, NULL, NULL, NULL};

fvm_nodal_t  *_nodal_src = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Return size associated with elements of a given restart type.
 *
 * parameters:
 *   val_type        <-- data type
 *
 * return:
 *   associated element type, in bytes
 *----------------------------------------------------------------------------*/

static size_t
_type_size(cs_restart_val_type_t   val_type)
{
  size_t  retval = 0;

  switch (val_type) {
  case CS_TYPE_char:
    retval = 1;
    break;
  case CS_TYPE_int:
    retval = sizeof(int);
    break;
  case CS_TYPE_cs_gnum_t:
    retval = sizeof(cs_gnum_t);
    break;
  case CS_TYPE_cs_real_t:
    retval = sizeof(cs_real_t);
    break;
  default:
    assert(0);
  }

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Use P0 interpolation (projection) from source to destination
 *        mesh entity.
 *
 * \param[in]       locator          associated locator
 * \param[in]       n_location_vals  number of values per location (interlaced)
 * \param[in]       val_type         value type
 * \param[in]       val_src          array of source values
 * \param[out]      val              array of values
 */
/*----------------------------------------------------------------------------*/

static void
_interpolate_p0(ple_locator_t          *locator,
                int                     n_location_vals,
                cs_restart_val_type_t   val_type,
                const void             *val_src,
                void                   *val)
{
  const unsigned char *_val_src = (const unsigned char *)val_src;

  size_t type_size = _type_size(val_type);
  size_t loc_size = type_size*n_location_vals;

  size_t  n_dist = ple_locator_get_n_dist_points(locator);
  const cs_lnum_t  *dist_loc = ple_locator_get_dist_locations(locator);

  /* Prepare and send data */

  unsigned char  *send_var;
  BFT_MALLOC(send_var, n_dist*loc_size, unsigned char);

  for (size_t i = 0; i < n_dist; i++) {
    const unsigned char  *src = _val_src + dist_loc[i]*loc_size;
    unsigned char  *dest = send_var + i*loc_size;
    for (size_t j = 0; j < loc_size; j++)
      dest[j] = src[j];
  }

  ple_locator_exchange_point_var_all(locator,
                                     send_var,
                                     val,
                                     NULL,
                                     type_size,
                                     n_location_vals,
                                     0);

  BFT_FREE(send_var);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Interolate at vertices from source to destination mesh entity.
 *
 * \param[in]       locator          associated locator
 * \param[in]       n_location_vals  number of values per location (interlaced)
 * \param[in]       val_type         value type
 * \param[in]       val_src          array of source values
 * \param[out]      val              array of values
 */
/*----------------------------------------------------------------------------*/

static void
_interpolate_vtx(ple_locator_t          *locator,
                 int                     n_location_vals,
                 cs_restart_val_type_t   val_type,
                 const void             *val_src,
                 void                   *val)
{
  cs_assert(val_type == CS_TYPE_cs_real_t);

  size_t type_size = _type_size(val_type);
  size_t loc_size = type_size*n_location_vals;

  size_t  n_dist = ple_locator_get_n_dist_points(locator);
  const cs_lnum_t  *dist_loc = ple_locator_get_dist_locations(locator);
  const ple_coord_t  *point_coords = ple_locator_get_dist_coords(locator);

  /* Prepare and send data */

  unsigned char  *send_var;
  BFT_MALLOC(send_var, n_dist*loc_size, unsigned char);

  cs_assert(_nodal_src != NULL);

  fvm_interpolate_vtx_data(_nodal_src,
                           3,                /* entity_dim */
                           n_location_vals,  /* data_dim */
                           n_dist,
                           dist_loc,
                           point_coords,
                           (const cs_real_t *)val_src,
                           (cs_real_t *)send_var);

  ple_locator_exchange_point_var_all(locator,
                                     send_var,
                                     val,
                                     NULL,
                                     type_size,
                                     n_location_vals,
                                     0);

  BFT_FREE(send_var);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read a section with interpolation.
 *
 * Note: if the context pointer is non-NULL, it must point to valid data
 * when the selection function is called, so that value or structure should
 * in general not be temporary.
 *
 * \param[in]       restart          associated restart file pointer
 * \param[in, out]  context          pointer to optional (untyped) value or
 *                                   structure.
 * \param[in]       sec_name         section name
 * \param[in]       location_id      id of corresponding location
 * \param[in]       n_location_vals  number of values per location (interlaced)
 * \param[in]       val_type         value type
 * \param[out]      val              array of values
 *
 * \return  0 (CS_RESTART_SUCCESS) in case of success,
 *          or error code (CS_RESTART_ERR_xxx) in case of error
 */
/*----------------------------------------------------------------------------*/

static int
_read_section_interpolate(cs_restart_t           *restart,
                          void                   *context,
                          const char             *sec_name,
                          int                     location_id,
                          int                     n_location_vals,
                          cs_restart_val_type_t   val_type,
                          void                   *val)
{
  CS_UNUSED(context);

  int retval = CS_RESTART_ERR_EXISTS;

  if (location_id == CS_MESH_LOCATION_NONE)
    retval = _read_section_f(restart,
                             context,
                             sec_name,
                             location_id,
                             n_location_vals,
                             val_type,
                             val);

  else {

    ple_locator_t  *locator = NULL;
    if (location_id < 5)
      locator = _locator[location_id - 1];

    if (locator == NULL)
      return CS_RESTART_ERR_NO_MAP;

    const cs_lnum_t n_src_elts
      = cs_restart_get_n_location_elts(restart, location_id);

    size_t type_size = _type_size(val_type);
    size_t loc_size = type_size*n_location_vals;

    unsigned char *read_buffer;
    BFT_MALLOC(read_buffer, n_src_elts*loc_size, unsigned char);

    retval = _read_section_f(restart,
                             context,
                             sec_name,
                             location_id,
                             n_location_vals,
                             val_type,
                             read_buffer);

    if (retval == CS_RESTART_SUCCESS) {

      if (location_id < CS_MESH_LOCATION_VERTICES)
        _interpolate_p0(locator,
                        n_location_vals,
                        val_type,
                        read_buffer,
                        val);
      else if (location_id == CS_MESH_LOCATION_VERTICES) {
        if (   _apply_mesh_deformation
            && strncmp(sec_name, "mesh_displacement::vals::", 25) == 0) {
          /* New mesh located relative to deformed mesh, so displacement
             reset to zero in this case. */
          cs_real_t v_0[] = {0, 0, 0};
          cs_array_real_set_vector(cs_glob_mesh->n_vertices,
                                   v_0,
                                   val);
        }
        else
          _interpolate_vtx(locator,
                           n_location_vals,
                           val_type,
                           read_buffer,
                           val);
      }

    }

    BFT_FREE(read_buffer);
  }

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read mesh deformation
 */
/*----------------------------------------------------------------------------*/

static void
_read_mesh_deformation(cs_mesh_t  *m)
{
  cs_restart_t *r
    = cs_restart_create("auxiliary.csc", "restart", CS_RESTART_MODE_READ);

  if (r == NULL)
    return;

  cs_real_3_t *v_disp;
  BFT_MALLOC(v_disp, m->n_vertices, cs_real_3_t);

  int retcode = cs_restart_read_section(r,
                                        "mesh_displacement::vals::0",
                                        CS_MESH_LOCATION_VERTICES,
                                        3,
                                        CS_TYPE_cs_real_t,
                                        v_disp);

  if (retcode == CS_RESTART_SUCCESS) {

    for (cs_lnum_t i = 0; i < m->n_vertices; i++) {
      for (cs_lnum_t j= 0; j < 3; j++)
        m->vtx_coord[i*3 + j] += v_disp[i][j];
    }

  }

  BFT_FREE(v_disp);

  cs_restart_destroy(&r);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Indicate restart files should be mapped to a given mesh input
 *
 * \param[in]  mesh_path           path to mesh input
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_map_set_mesh_input(const char  *mesh_path)
{
  size_t n = strlen(mesh_path);
  BFT_REALLOC(_mesh_input_path, n + 1, char);

  strncpy(_mesh_input_path, mesh_path, n + 1);
  _mesh_input_path[n] = '\0';
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set options relative to restart file mapping to a given mesh input.
 *
 * \param[in]  apply_mesh_deformation  apply mesh deformation from upstream
 *                                     computation (if present) so as to map
 *                                     to final, and not initial mesh shape.
 * \param[in]  tolerance_base          associated base tolerance (used for
 *                                     bounding box check only, not for
 *                                     location test)
 * \param[in]  tolerance_fraction      associated fraction of element bounding
 *                                     boxes added to tolerance
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_map_set_options(bool   apply_mesh_deformation,
                           float  tolerance_base,
                           float  tolerance_fraction)
{
  _apply_mesh_deformation = apply_mesh_deformation;

  _tolerance[0] = tolerance_base;
  _tolerance[1] = tolerance_fraction;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Indicate whether location for restart file mapping is needed at
 *         cells or vertices.
 *
 * By default, mapping is done for cell-based quantities, but not for
 * vertex-based quantities.
 *
 * Mapping of quantities at faces or particles is not handled yet, but will
 * use the cell-center or vertex based mappings in the future in all cases:
 * - interior faces may not be aligned with previous faces, so some sort of
 *   cell-based interpolation will be required
 * - boundary faces can use the boundary face / cell adjacency to avoid an
 *   additional mapping
 * - for particles, as the previous location is stored based on cell ids,
 *   updating particle locations will require locating them in the matching
 *   cell and completing a trajectory adjustment (direct location should be
 *   avoided in case of concave boundaries).
 *
 * \param[in]  map_cell_centers    locate cell centers in the previous mesh.
 * \param[in]  map_vertices        locate vertices in the previous mesh.
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_map_set_locations(bool map_cell_centers,
                             bool map_vertices)
{
  _need_locator[0] = map_cell_centers;
  _need_locator[1] = false;
  _need_locator[2] = false;
  _need_locator[3] = map_vertices;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build mapping of restart files to different mesh if defined.
 *
 * \param[in]  need_vertices  indicate if location at vertices is needed
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_map_build(void)
{
  if (_mesh_input_path == NULL)
    return;

  int t_restart_id = cs_timer_stats_id_by_name("checkpoint_restart_stage");
  int t_top_id = cs_timer_stats_switch(t_restart_id);

  /* Stash (protect) mesh to read older mesh; should not be necessary
     for reading mesh, but required for older restart, and
     may be safer at this stage */

  cs_mesh_t *m_c = cs_glob_mesh;
  cs_glob_mesh = NULL;

  /* Read previous mesh */

  fvm_nodal_t  *nm = NULL;

  {
    /* Read mesh */

    cs_mesh_t *m = cs_mesh_create();

    cs_mesh_builder_t *mb_s = cs_glob_mesh_builder;
    cs_glob_mesh_builder = NULL;

    cs_mesh_builder_t *mb = cs_mesh_builder_create();

    cs_preprocessor_data_add_file(_mesh_input_path, 0, NULL, NULL);
    /* As currently we restart from an existing file, we set the
     * "ignore_cartesian" option to true. */
    cs_preprocessor_data_read_headers(m, mb, true);
    cs_preprocessor_data_read_mesh(m, mb, true);

    m->n_b_faces_all = m->n_b_faces;
    m->n_g_b_faces_all = m->n_g_b_faces;

    cs_mesh_builder_destroy(&mb);
    cs_glob_mesh_builder = mb_s;

    /* Set reference numbering for restart */

    cs_restart_add_location_ref("cells",
                                m->n_g_cells, m->n_cells,
                                m->global_cell_num);
    cs_restart_add_location_ref("interior_faces",
                                m->n_g_i_faces, m->n_i_faces,
                                m->global_i_face_num);
    cs_restart_add_location_ref("boundary_faces",
                                m->n_g_b_faces, m->n_b_faces,
                                m->global_b_face_num);
    cs_restart_add_location_ref("vertices",
                                m->n_g_vertices, m->n_vertices,
                                m->global_vtx_num);

    cs_glob_mesh = m;

    if (_apply_mesh_deformation)
      _read_mesh_deformation(m);

    /* Build FVM mesh from previous mesh */

    nm = cs_mesh_connect_cells_to_nodal(m,
                                        "restart_mesh",
                                        false,
                                        m->n_cells,
                                        NULL);

    fvm_nodal_make_vertices_private(nm);

    /* Destroy temporary mesh structures */

    m = cs_mesh_destroy(m);
  }

  cs_glob_mesh = m_c;

  /* Now build locator */
  /*-------------------*/

  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  int options[PLE_LOCATOR_N_OPTIONS];
  for (int i = 0; i < PLE_LOCATOR_N_OPTIONS; i++)
    options[i] = 0;
  options[PLE_LOCATOR_NUMBERING] = 0; /* base 0 numbering */

  cs_lnum_t n_points[4] = {m_c->n_cells, 0, 0, m_c->n_vertices};
  const cs_real_t *point_coords[4] = {mq->cell_cen, NULL, NULL, m_c->vtx_coord};

  for (int i = 0; i < 4; i++) {

    if (_need_locator[i] == false) continue;

#if defined(PLE_HAVE_MPI)
    _locator[i] = ple_locator_create(cs_glob_mpi_comm,
                                     cs_glob_n_ranks,
                                     0);
#else
    _locator[i] = ple_locator_create();
#endif

    /* For cells (or faces), locate relative to parent mesh.
       For vertices, locate relative to sections, as this will be better
       adapted to the fvm_point_location_interpolate_vtx_data call. */

    ple_mesh_elements_locate_t  *mesh_elements_locate_f
      = cs_coupling_point_in_mesh_p;
    if (i == 3)
      mesh_elements_locate_f = cs_coupling_point_in_mesh;

    ple_locator_set_mesh(_locator[i],
                         nm,
                         options,
                         _tolerance[0],
                         _tolerance[1],
                         3, /* dim */
                         n_points[i],
                         NULL,
                         NULL, /* point_tag */
                         point_coords[i],
                         NULL, /* distance */
                         cs_coupling_mesh_extents,
                         mesh_elements_locate_f);

    /* Shift from 1-base to 0-based locations */

    ple_locator_shift_locations(_locator[i], -1);

  }

  /* Nodal mesh may not be needed anymore */

  if (_need_locator[3]) {
    _nodal_src = nm;
    nm = NULL;
  }
  else
    nm = fvm_nodal_destroy(nm);

  /* Set associated read function if not already set */

  if (_read_section_f == NULL) {
    _read_section_f
      = cs_restart_set_read_section_func(_read_section_interpolate);
  }

  cs_timer_stats_switch(t_top_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free restart file mapping to different mesh.
 *
 * Revert restart reading to default behavior.
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_map_free(void)
{
  BFT_FREE(_mesh_input_path);
  _tolerance[0] = 0;
  _tolerance[1] = 0.1;

  if (_nodal_src != NULL)
    _nodal_src = fvm_nodal_destroy(_nodal_src);

  if (_read_section_f != NULL) {
    (void)cs_restart_set_read_section_func(_read_section_f);
    _read_section_f = NULL;

    cs_restart_clear_locations_ref();
  }

  double loc_times[4] = {0, 0, 0, 0};

  for (int i = 0; i < 4; i++ ) {

    if (_locator[i] != NULL) {

      double _loc_times[4];

      ple_locator_get_times(_locator[i],
                            _loc_times,
                            NULL,
                            _loc_times+1,
                            NULL);
      ple_locator_get_comm_times(_locator[i],
                                 _loc_times+2,
                                 NULL,
                                 _loc_times+3,
                                 NULL);

      for (int j = 0; j < 4; j++)
        loc_times[j] += _loc_times[j];

    }

  }

  if (cs_glob_n_ranks < 2)
    cs_log_printf
      (CS_LOG_PERFORMANCE,
       _("\n"
         "Restart mapping\n"
         "                                 "
         "  location time:                 %12.3f\n"
         "    communication and wait:      %12.3f\n"
         "  variable exchange time:        %12.3f\n"
         "    communication and wait:      %12.3f\n\n"),
       loc_times[0], loc_times[2], loc_times[1], loc_times[3]);

  else {

    double max_times[4], min_times[4], mean_times[4];

    for (int i = 0; i < 4; i++) {
      max_times[i] = loc_times[i];
      min_times[i] = loc_times[i];
      mean_times[i] = loc_times[i];
    }
    cs_parall_min(4, CS_DOUBLE, min_times);
    cs_parall_max(4, CS_DOUBLE, max_times);
    cs_parall_sum(4, CS_DOUBLE, mean_times);
    for (int i = 0; i < 4; i++)
      mean_times[i] /= cs_glob_n_ranks;

    cs_log_printf
      (CS_LOG_PERFORMANCE,
       _("\n"
         "Restart mapping\n"
         "                                 "
         "        mean      minimum     maximum\n"
         "  location time:                 %12.3f %12.3f %12.3f\n"
         "    communication and wait:      %12.3f %12.3f %12.3f\n"
         "  variable exchange time:        %12.3f %12.3f %12.3f\n"
         "    communication and wait:      %12.3f %12.3f %12.3f\n\n"),
       mean_times[0], min_times[0], max_times[0],
       mean_times[2], min_times[2], max_times[2],
       mean_times[1], min_times[1], max_times[1],
       mean_times[3], min_times[3], max_times[3]);

  }

  cs_log_separator(CS_LOG_PERFORMANCE);

  for (int i = 0; i < 4; i++ ) {
    if (_locator[i] != NULL)
      _locator[i] = ple_locator_destroy(_locator[i]);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
