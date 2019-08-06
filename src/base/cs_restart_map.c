/*============================================================================
 * Checkpoint / restart extension to mapped meshes
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"

#include "fvm_nodal.h"
#include "ple_locator.h"

#include "cs_io.h"
#include "cs_coupling.h"
#include "cs_mesh.h"
#include "cs_mesh_connect.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_preprocessor_data.h"
#include "cs_timer.h"
#include "cs_timer_stats.h"

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

/* Variables to be saved between definition of mapping and
   build of tht mapping */

static char *_mesh_input_path = NULL;
static float _tolerance[2] = {0, 0.1};

/* Save previous section read function */
static cs_restart_read_section_t   *_read_section_f  = NULL;

static  ple_locator_t  *_locator = NULL;  /* PLE locator for restart */

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
  case CS_TYPE_cs_int_t:
    retval = sizeof(cs_int_t);
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

  ple_locator_exchange_point_var(_locator,
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

  else if (location_id == CS_MESH_LOCATION_CELLS) {

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

    if (retval == CS_RESTART_SUCCESS)
      _interpolate_p0(_locator,
                      n_location_vals,
                      val_type,
                      read_buffer,
                      val);

    BFT_FREE(read_buffer);
  }

  return retval;
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Indicate restart files should be mapped to a given mesh input
 *
 * \param[in]  mesh_path           path to mesh input
 * \param[in]  tolerance_base      associated base tolerance (used for bounding
 *                                 box check only, not for location test)
 * \param[in]  tolerance_fraction  associated fraction of element bounding
 *                                 boxes added to tolerance
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_map_set_mesh_input(const char  *mesh_path)
{
  size_t n = strlen(mesh_path);
  BFT_REALLOC(_mesh_input_path, n + 1, char);

  strncpy(_mesh_input_path, mesh_path, n);
  _mesh_input_path[n] = '\0';
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set options relative to restart file mapping to a given mesh input.
 *
 * \param[in]  tolerance_base      associated base tolerance (used for bounding
 *                                 box check only, not for location test)
 * \param[in]  tolerance_fraction  associated fraction of element bounding
 *                                 boxes added to tolerance
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_map_set_options(float  tolerance_base,
                           float  tolerance_fraction)
{
  _tolerance[0] = tolerance_base;
  _tolerance[1] = tolerance_fraction;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build mapping of restart files to different mesh if defined.
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
    cs_preprocessor_data_read_headers(m, mb);
    cs_preprocessor_data_read_mesh(m, mb);

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

    /* Build FVM mesh from previous mesh */

    nm = cs_mesh_connect_cells_to_nodal(m,
                                        "restart_mesh",
                                        false,
                                        m->n_cells,
                                        NULL);

    fvm_nodal_make_vertices_private(nm);

    /* Destroy temporary mesh structures */

    cs_glob_mesh = m;
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

#if defined(PLE_HAVE_MPI)
  _locator = ple_locator_create(cs_glob_mpi_comm,
                                cs_glob_n_ranks,
                                0);
#else
  _locator = ple_locator_create();
#endif

  ple_locator_set_mesh(_locator,
                       nm,
                       options,
                       _tolerance[0],
                       _tolerance[1],
                       3, /* dim */
                       m_c->n_cells,
                       NULL,
                       NULL, /* point_tag */
                       mq->cell_cen,
                       NULL, /* distance */
                       cs_coupling_mesh_extents,
                       cs_coupling_point_in_mesh_p);

  /* Shift from 1-base to 0-based locations */

  ple_locator_shift_locations(_locator, -1);

  /* Nodal mesh is not needed anymore */

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

  if (_read_section_f != NULL) {
    (void)cs_restart_set_read_section_func(_read_section_f);
    _read_section_f = NULL;

    cs_restart_clear_locations_ref();
  }

  _locator = ple_locator_destroy(_locator);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
