/*============================================================================
 * Define postprocessing output.
 *============================================================================*/

/* VERS */

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"

#include "fvm_writer.h"

#include "cs_base.h"
#include "cs_field.h"
#include "cs_mesh.h"
#include "cs_selector.h"
#include "cs_probe.h"
#include "cs_post.h"
#include "cs_time_plot.h"
#include "cs_order.h"
#include "cs_part_to_block.h"
#include "fvm_nodal_append.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local (user defined) function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Write space-filling curves for main mesh
 *
 * parameters:
 *   writer <-- FVM writer
 *---------------------------------------------------------------------------*/

static void
_cs_post_write_sfc_serial(fvm_writer_t   *writer)
{
  fvm_io_num_sfc_t  sfc_id;
  cs_lnum_t  i, j, k;

  cs_lnum_t *connect = NULL, *order = NULL;
  double *coords = NULL, *val = NULL;
  fvm_nodal_t *nm = NULL;
  fvm_io_num_t *io_num = NULL;

  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_lnum_t n_edges = m->n_cells - 1;
  const cs_gnum_t *cell_gnum = NULL;
  const double *var_ptr[1] = {NULL};

  BFT_MALLOC(order, m->n_cells, cs_lnum_t);
  BFT_MALLOC(val, m->n_cells, double);
  BFT_MALLOC(coords, m->n_cells*3, double);

  /* Loop on space-filling curve types */

  for (sfc_id = FVM_IO_NUM_SFC_MORTON_BOX;
       sfc_id <= FVM_IO_NUM_SFC_HILBERT_CUBE;
       sfc_id++) {

    BFT_MALLOC(connect, n_edges*2, cs_lnum_t);

    io_num = fvm_io_num_create_from_sfc(mq->cell_cen,
                                        3,
                                        m->n_cells,
                                        sfc_id);

    cell_gnum = fvm_io_num_get_global_num(io_num);

    cs_order_gnum_allocated(NULL, cell_gnum, order, m->n_cells);

    for (i = 0; i < m->n_cells; i++) {
      j = order[i];
      for (k = 0; k < 3; k++)
        coords[i*3 + k] = mq->cell_cen[j*3 + k];
      val[i] = i+1;
    }

    for (i = 0; i < n_edges; i++) {
      connect[i*2] = i+1;
      connect[i*2+1] = i+2;
    }

    cell_gnum = NULL;
    fvm_io_num_destroy(io_num);

    nm = fvm_nodal_create(fvm_io_num_sfc_type_name[sfc_id], 3);

    fvm_nodal_append_by_transfer(nm,
                                 m->n_cells - 1,
                                 FVM_EDGE,
                                 NULL,
                                 NULL,
                                 NULL,
                                 connect,
                                 NULL);

    fvm_nodal_set_shared_vertices(nm, coords);

    fvm_writer_export_nodal(writer, nm);

    var_ptr[0] = val;

    fvm_writer_export_field(writer,
                            nm,
                            _("order"),
                            FVM_WRITER_PER_NODE,
                            1,
                            CS_INTERLACE,
                            0,
                            0,
                            CS_DOUBLE,
                            -1,
                            0.0,
                            (const void * *)var_ptr);

    fvm_nodal_destroy(nm);
  }

  /* Free memory */

  BFT_FREE(order);
  BFT_FREE(val);
  BFT_FREE(coords);
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Write space-filling curves for main mesh
 *
 * parameters:
 *   writer <-- FVM writer
 *---------------------------------------------------------------------------*/

static void
_cs_post_write_sfc_parall(fvm_writer_t  *writer)
{
  fvm_io_num_sfc_t  sfc_id;
  cs_lnum_t  i;

  cs_lnum_t *connect = NULL, *order = NULL;
  cs_gnum_t *vtx_gnum = NULL, *edge_gnum = NULL;
  double *val = NULL;
  cs_coord_t *coords = NULL;
  fvm_nodal_t *nm = NULL;
  fvm_io_num_t *io_num = NULL;

  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_gnum_t *cell_gnum = NULL;
  const double  *var_ptr[1] = {NULL};

  /* Loop on space-filling curve types */

  for (sfc_id = FVM_IO_NUM_SFC_MORTON_BOX;
       sfc_id <= FVM_IO_NUM_SFC_HILBERT_CUBE;
       sfc_id++) {

    cs_lnum_t block_size = 0;
    cs_lnum_t n_edges = 0;
    cs_block_dist_info_t bi;
    cs_part_to_block_t *d = NULL;

    io_num = fvm_io_num_create_from_sfc(mq->cell_cen,
                                        3,
                                        m->n_cells,
                                        sfc_id);

    cell_gnum = fvm_io_num_get_global_num(io_num);

    /* Distribute to blocks so that edge connectivity is trivial */

    bi = cs_block_dist_compute_sizes(cs_glob_rank_id,
                                     cs_glob_n_ranks,
                                     0,
                                     0,
                                     m->n_g_cells);

    d = cs_part_to_block_create_by_gnum(cs_glob_mpi_comm,
                                        bi,
                                        m->n_cells,
                                        cell_gnum);

    block_size = (bi.gnum_range[1] - bi.gnum_range[0]);

    if (block_size > 0) {
      BFT_MALLOC(connect, block_size*2, cs_lnum_t);
      BFT_MALLOC(val, block_size+1, double);
      BFT_MALLOC(coords, (block_size+1)*3, double);
      BFT_MALLOC(vtx_gnum, block_size+1, cs_gnum_t);
      BFT_MALLOC(edge_gnum, block_size, cs_gnum_t);
    }

    /* Distribute blocks on ranks */

    cs_part_to_block_copy_array(d,
                                CS_DOUBLE,
                                3,
                                mq->cell_cen,
                                coords);

    cell_gnum = NULL;
    fvm_io_num_destroy(io_num);

    cs_part_to_block_destroy(&d);

    /* Add vertex for connectivity with next rank */

    if (block_size > 0) {
      MPI_Status status;
      int prev_rank = cs_glob_rank_id - bi.rank_step;
      int next_rank = cs_glob_rank_id + bi.rank_step;
      if (prev_rank < 0)
        prev_rank = MPI_PROC_NULL;
      if (bi.gnum_range[1] > m->n_g_cells)
        next_rank = MPI_PROC_NULL;
      MPI_Sendrecv(coords, 3, MPI_DOUBLE, prev_rank, 0,
                   coords + 3*block_size, 3, MPI_DOUBLE, next_rank, 0,
                   cs_glob_mpi_comm, &status);
    }

    for (i = 0; i < block_size; i++) {
      vtx_gnum[i] = bi.gnum_range[0] + i;
      if (vtx_gnum[i] < m->n_g_cells) {
        connect[n_edges*2] = i+1;
        connect[n_edges*2+1] = i+2;
        edge_gnum[n_edges] = vtx_gnum[i];
        n_edges++;
      }
      val[i] = vtx_gnum[i];
    }
    if (block_size > 0) {
      vtx_gnum[block_size] = bi.gnum_range[0] + block_size;
      val[block_size] = vtx_gnum[block_size];
    }

    BFT_FREE(order);

    nm = fvm_nodal_create(fvm_io_num_sfc_type_name[sfc_id], 3);

    fvm_nodal_append_by_transfer(nm,
                                 n_edges,
                                 FVM_EDGE,
                                 NULL,
                                 NULL,
                                 NULL,
                                 connect,
                                 NULL);

    connect = NULL;

    fvm_nodal_set_shared_vertices(nm, coords);

    fvm_nodal_init_io_num(nm, edge_gnum, 1);
    fvm_nodal_init_io_num(nm, vtx_gnum, 0);

    fvm_writer_export_nodal(writer, nm);

    var_ptr[0] = val;

    fvm_writer_export_field(writer,
                            nm,
                            _("order"),
                            FVM_WRITER_PER_NODE,
                            1,
                            CS_INTERLACE,
                            0,
                            0,
                            CS_DOUBLE,
                            -1,
                            0.0,
                            (const void * *)var_ptr);

    /* Free memory */

    if (block_size > 0) {
      BFT_FREE(val);
      BFT_FREE(coords);
      BFT_FREE(vtx_gnum);
      BFT_FREE(edge_gnum);
    }

    fvm_nodal_destroy(nm);

  }
}

#endif /* defined(HAVE_MPI) */

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define post-processing writers.
 *
 * The default output format and frequency may be configured, and additional
 * post-processing writers allowing outputs in different formats or with
 * different format options and output frequency than the main writer may
 * be defined.
 *----------------------------------------------------------------------------*/

void
cs_user_postprocess_writers(void)
{
  /* Every writer has a strictly positive or negative id. Negative ids
   * are for predefined writers, positive ids for user writers.
   * All predefined writers use the settings from writer -1, and
   * redefining that writer here allows changing from the default or GUI
   * settings.
   *
   * Defining or configuring a writer is done by calling the
   * cs_post_define_writer() function, whose arguments are:
   *   writer_id     <-- number of writer to create (< 0 reserved, > 0 for user)
   *   case_name     <-- associated case name
   *   dir_name      <-- associated directory name
   *   fmt_name      <-- associated format name
   *   fmt_opts      <-- associated format options string
   *   time_dep      <-- FVM_WRITER_FIXED_MESH if mesh definitions are fixed,
   *                     FVM_WRITER_TRANSIENT_COORDS if coordinates change,
   *                     FVM_WRITER_TRANSIENT_CONNECT if connectivity changes
   *   output_at_end <-- force output at calculation end if not 0
   *   frequency_n   <-- default output frequency in time-steps, or < 0
   *   frequency_t   <-- default output frequency in seconds, or < 0
   *                     (has priority over frequency_n)
   *
   * Allowed output format names: "EnSight Gold", "MED", or "CGNS".
   * (EnSight output is built-in; MED or CGNS are only available if the
   * code was built with these optional libraries)
   *
   * An output options string may contain options (separated by whitespace
   * or commas) from the following list:
   *   'text'              (text format, for EnSight)
   *   'big_endian'        (forces binary EnSight output to 'big-endian' mode)
   *   'adf'               (use ADF file type, for CGNS)
   *   'hdf5'              (force HDF5 file type, usually the default for CGNS)
   *   'discard_polygons'  (ignore polygon-type faces)
   *   'discard_polyhedra' (ignore polyhedron-type cells)
   *   'divide_polygons'   (subdivides polygon-type faces)
   *   'divide_polyhedra'  (subdivides polyhedron-type cells)
   *   'separate_meshes'   (use separate outputs for different meshes) */

  /* Set time plot file writer flush behavior defaults.
   *   flush_wtime     <-- elapsed time interval between file flushes;
   *                       if < 0, no forced flush
   *   n_buffer_steps  <-- number of time steps in output buffer if
   *                       file is not to be kept open */
}

/*----------------------------------------------------------------------------
 * Define post-processing meshes.
 *
 * The main post-processing meshes may be configured, and additional
 * post-processing meshes may be defined as a subset of the main mesh's
 * cells or faces (both interior and boundary).
 *----------------------------------------------------------------------------*/

void
cs_user_postprocess_meshes(void)
{
  /* In this example, a specific, temporary writer is built,
     using the defaults (based on writer -1), and
     edge meshes illustrating the various space-filling curve
     possibilities are output using this writer. */

  fvm_writer_t *w = NULL;

  /* Create default writer */

  w = fvm_writer_init("SFC",
                      "postprocessing",
                      cs_post_get_default_format(),
                      cs_post_get_default_format_options(),
                      FVM_WRITER_FIXED_MESH);

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1)
    _cs_post_write_sfc_parall(w);
#endif

  if (cs_glob_n_ranks == 1)
    _cs_post_write_sfc_serial(w);

  fvm_writer_finalize(w);
}

/*----------------------------------------------------------------------------
 * Override default frequency or calculation end based output.
 *
 * This allows fine-grained control of activation or deactivation,
 *
 * parameters:
 *   nt_max_abs <-- maximum time step number
 *   nt_cur_abs <-- current time step number
 *   t_cur_abs  <-- absolute time at the current time step
 *----------------------------------------------------------------------------*/

void
cs_user_postprocess_activate(int     nt_max_abs,
                             int     nt_cur_abs,
                             double  t_cur_abs)
{
  /* Use the cs_post_activate_writer() function to force the
   * "active" or "inactive" flag for a specific writer or for all
   * writers for the current time step.

   * the parameters for cs_post_activate_writer() are:
   *   writer_id <-- writer id, or 0 for all writers
   *   activate  <-- false to deactivate, true to activate */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
