/*============================================================================
 * Manage the exchange of data between Code_Saturne and the pre-processor
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2012 EDF S.A.

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
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_error.h>
#include <bft_mem.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_io_num.h"
#include "fvm_periodicity.h"

#include "cs_base.h"
#include "cs_block_dist.h"
#include "cs_block_to_part.h"
#include "cs_file.h"
#include "cs_interface.h"
#include "cs_mesh.h"
#include "cs_io.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_preprocessor_data.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/* Structure used for building mesh structure */
/* ------------------------------------------ */

typedef struct {

  const char         *filename;   /* File name */
  cs_file_off_t       offset;     /* File offsets for re-opening */
  const double       *matrix;     /* Coordinate transformation matrix */

  size_t              n_group_renames;
  const char  *const *old_group_names;
  const char  *const *new_group_names;

  /* Single allocation for all data */

  size_t              data_size;
  unsigned char      *data;

} _mesh_file_info_t;

/* Structure used for building mesh structure */
/* ------------------------------------------ */

typedef struct {

  /* File info */

  int                 n_files;
  _mesh_file_info_t  *file_info;

  /* Face-related dimensions */

  cs_gnum_t   n_g_faces;
  cs_gnum_t   n_g_face_connect_size;

  /* Temporary dimensions necessary for multiple inputs */

  int         *gc_id_shift;

  int          n_perio_read;
  cs_lnum_t    n_cells_read;
  cs_lnum_t    n_faces_read;
  cs_lnum_t    n_faces_connect_read;
  cs_lnum_t    n_vertices_read;

  cs_gnum_t    n_g_cells_read;
  cs_gnum_t    n_g_faces_read;
  cs_gnum_t    n_g_faces_connect_read;
  cs_gnum_t    n_g_vertices_read;

  /* Temporary mesh data */

  int           read_cell_rank;
  int          *cell_rank;

  cs_gnum_t    *face_cells;
  cs_lnum_t    *face_vertices_idx;
  cs_gnum_t    *face_vertices;
  cs_int_t     *cell_gc_id;
  cs_int_t     *face_gc_id;
  cs_real_t    *vertex_coords;

  /* Periodic features */

  int           n_perio;               /* Number of periodicities */
  int          *periodicity_num;       /* Periodicity numbers */
  cs_lnum_t    *n_per_face_couples;    /* Nb. face couples per periodicity */
  cs_gnum_t    *n_g_per_face_couples;  /* Global nb. couples per periodicity */

  cs_gnum_t   **per_face_couples;      /* Periodic face couples list. */

  /* Block ranges for parallel distribution */

  cs_block_dist_info_t   cell_bi;     /* Block info for cell data */
  cs_block_dist_info_t   face_bi;     /* Block info for face data */
  cs_block_dist_info_t   vertex_bi;   /* Block info for vertex data */
  cs_block_dist_info_t  *per_face_bi; /* Block info for parallel face couples */

} _mesh_reader_t;

typedef double  _vtx_coords_t[3];

/*============================================================================
 *  Global variables
 *============================================================================*/

static bool             _use_sfc = true;
static fvm_io_num_sfc_t _sfc_type = FVM_IO_NUM_SFC_MORTON_BOX;

static _mesh_reader_t *_cs_glob_mesh_reader = NULL;

#if defined(WIN32) || defined(_WIN32)
static const char _dir_separator = '\\';
#else
static const char _dir_separator = '/';
#endif

/* Definitions of file to read */

int _n_mesh_files = 0;
int _n_max_mesh_files = 0;
_mesh_file_info_t  *_mesh_file_info = NULL;

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Return minimum size required to align data with the base pointer size.
 *
 * parameters:
 *   min_size <-- minimum data size
 *
 * returns:
 *   the data size requirde including alignment
 *----------------------------------------------------------------------------*/

static inline size_t
_align_size(size_t  min_size)
{
  const size_t align = (sizeof(void *))/sizeof(unsigned char);
  return (min_size + (align-1) - ((min_size - 1) % align));
}

/*----------------------------------------------------------------------------
 * Define defaul input data in nothing has been specified by the user.
 *----------------------------------------------------------------------------*/

static void
_set_default_input_if_needed(void)
{
  const char input_default[] = "mesh_input";

  if (_n_mesh_files == 0) {

    if (cs_file_isreg(input_default))
      cs_preprocessor_data_add_file(input_default, 0, NULL, NULL);

    else if (cs_file_isdir(input_default)) {
      int i;
      char **dir_files = cs_file_listdir(input_default);
      for (i = 0; dir_files[i] != NULL; i++) {
        char *tmp_name = NULL;
        BFT_MALLOC(tmp_name,
                   strlen(input_default) + 1 + strlen(dir_files[i]) + 1,
                   char);
        sprintf(tmp_name, "%s%c%s",
                input_default, _dir_separator, dir_files[i]);
        if (cs_file_isreg(tmp_name))
          cs_preprocessor_data_add_file(tmp_name, 0, NULL, NULL);
        BFT_FREE(tmp_name);
        BFT_FREE(dir_files[i]);
      }
      BFT_FREE(dir_files);
    }

    else
      bft_error(__FILE__, __LINE__, 0,
                _("No \"%s\" file or directory found."), input_default);
  }
}

/*----------------------------------------------------------------------------
 * Create an empty mesh reader helper structure.
 *
 * Property of the matching mesh file information structure is transferred
 * to the mesh reader.
 *
 * parameters:
 *   n_mesh_files   <-> number of associated mesh files
 *   mesh_file_info <-> array of mesh file information structures
 *
 * returns:
 *   A pointer to a mesh reader helper structure
 *----------------------------------------------------------------------------*/

static _mesh_reader_t *
_mesh_reader_create(int                 *n_mesh_files,
                    _mesh_file_info_t  **mesh_file_info)
{
  int i;
  _mesh_reader_t  *mr = NULL;

  BFT_MALLOC(mr, 1, _mesh_reader_t);

  memset(mr, 0, sizeof(_mesh_reader_t));

  /* Transfer ownership of mesh file info */

  mr->n_files = *n_mesh_files;
  mr->file_info = *mesh_file_info;

  BFT_REALLOC(mr->file_info, mr->n_files, _mesh_file_info_t);

  /* Setup remaining structure fields */

  *n_mesh_files = 0;
  *mesh_file_info = NULL;

  BFT_MALLOC(mr->gc_id_shift, mr->n_files, int);
  for (i = 0; i < mr->n_files; i++)
    mr->gc_id_shift[i] = 0;

  mr->n_g_faces = 0;
  mr->n_g_face_connect_size = 0;

  mr->n_perio_read = 0;
  mr->n_cells_read = 0;
  mr->n_faces_read = 0;
  mr->n_faces_connect_read = 0;
  mr->n_vertices_read = 0;

  mr->n_g_cells_read = 0;
  mr->n_g_faces_read = 0;
  mr->n_g_faces_connect_read = 0;

  mr->read_cell_rank = 0;

  mr->cell_rank = NULL;
  mr->face_cells = NULL;
  mr->face_vertices_idx = NULL;
  mr->face_vertices = NULL;
  mr->cell_gc_id = NULL;
  mr->face_gc_id = NULL;
  mr->vertex_coords = NULL;

  mr->n_perio = 0;
  mr->periodicity_num = NULL;
  mr->n_per_face_couples = NULL;
  mr->n_g_per_face_couples = NULL;
  mr->per_face_couples = NULL;

  mr->per_face_bi = NULL;

  return mr;
}

/*----------------------------------------------------------------------------
 * Destroy a mesh reader helper structure
 *
 * parameters:
 *   mr <-> pointer to a mesh reader helper
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

static void
_mesh_reader_destroy(_mesh_reader_t  **mr)
{
  int i;
  _mesh_reader_t *_mr = *mr;

  for (i = 0; i < _mr->n_files; i++) {
    _mesh_file_info_t  *f = _mr->file_info + i;
    BFT_FREE(f->data);
  }
  BFT_FREE(_mr->file_info);

  BFT_FREE(_mr->gc_id_shift);

  BFT_FREE(_mr->face_cells);
  BFT_FREE(_mr->face_vertices_idx);
  BFT_FREE(_mr->face_vertices);
  BFT_FREE(_mr->cell_gc_id);
  BFT_FREE(_mr->face_gc_id);
  BFT_FREE(_mr->vertex_coords);

  if (_mr->n_perio > 0) {
    if (_mr->per_face_couples != NULL) {
      for (i = 0; i < _mr->n_perio; i++)
        BFT_FREE(_mr->per_face_couples[i]);
    }
    BFT_FREE(_mr->per_face_couples);
    BFT_FREE(_mr->n_g_per_face_couples);
    BFT_FREE(_mr->n_per_face_couples);
    BFT_FREE(_mr->periodicity_num);
    BFT_FREE(_mr->per_face_bi);
  }

  BFT_FREE(*mr);
}

/*----------------------------------------------------------------------------
 * Add a periodicity to mesh->periodicities (fvm_periodicity_t *) structure.
 *
 * parameters:
 *   mesh       <-> mesh
 *   perio_type <-- periodicity type
 *   perio_num  <-- periodicity number (identifier)
 *   matrix     <-- transformation matrix using homogeneous coordinates
 *----------------------------------------------------------------------------*/

static void
_add_periodicity(cs_mesh_t *mesh,
                 cs_int_t   perio_type,
                 cs_int_t   perio_num,
                 cs_real_t  matrix[3][4])
{
  cs_int_t  i, j;
  double  _matrix[3][4];

  fvm_periodicity_type_t _perio_type = (fvm_periodicity_type_t)perio_type;

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 4; j++)
      _matrix[i][j] = matrix[i][j];
  }

  fvm_periodicity_add_by_matrix(mesh->periodicity,
                                perio_num,
                                _perio_type,
                                _matrix);
}

/*----------------------------------------------------------------------------
 * Set block ranges for parallel reads
 *
 * parameters:
 *   mesh <-- pointer to mesh structure
 *   mr   <-> mesh reader helper
 *----------------------------------------------------------------------------*/

static void
_set_block_ranges(cs_mesh_t       *mesh,
                  _mesh_reader_t  *mr)
{
  int i;

  int rank_id = cs_glob_rank_id;
  int n_ranks = cs_glob_n_ranks;

  /* Always build per_face_range in case of periodicity */

  if (mr->n_perio > 0) {
    BFT_REALLOC(mr->per_face_bi, mr->n_perio, cs_block_dist_info_t);
    memset(mr->per_face_bi, 0, sizeof(cs_block_dist_info_t)*mr->n_perio);
  }

  /* Set block sizes and ranges (useful for parallel mode) */

  mr->cell_bi = cs_block_dist_compute_sizes(rank_id,
                                            n_ranks,
                                            0,
                                            0,
                                            mesh->n_g_cells);

  mr->face_bi = cs_block_dist_compute_sizes(rank_id,
                                            n_ranks,
                                            0,
                                            0,
                                            mr->n_g_faces);

  mr->vertex_bi = cs_block_dist_compute_sizes(rank_id,
                                              n_ranks,
                                              0,
                                              0,
                                               mesh->n_g_vertices);

  for (i = 0; i < mr->n_perio; i++)
    mr->per_face_bi[i]
      = cs_block_dist_compute_sizes(rank_id,
                                    n_ranks,
                                    0,
                                    0,
                                    mr->n_g_per_face_couples[i]);
}

/*----------------------------------------------------------------------------
 * Read cell rank if available
 *
 * parameters:
 *   mesh <-- pointer to mesh structure
 *   mr   <-> mesh reader helper
 *   echo <-- echo (verbosity) level
 *----------------------------------------------------------------------------*/

static void
_read_cell_rank(cs_mesh_t       *mesh,
                _mesh_reader_t  *mr,
                long             echo)
{
  char file_name[64]; /* more than enough for
                         "partition/domain_number_<n_ranks>" */
  size_t  i;
  cs_io_sec_header_t  header;

  cs_io_t  *rank_pp_in = NULL;
  cs_lnum_t   n_ranks = 0;
  cs_gnum_t   n_elts = 0;
  cs_gnum_t   n_g_cells = 0;

  const char  *unexpected_msg = N_("Section of type <%s> on <%s>\n"
                                   "unexpected or of incorrect size");

  if (n_ranks == 1)
    return;

#if (__STDC_VERSION__ < 199901L)
  sprintf(file_name,
          "partition%cdomain_number_%d",
          _dir_separator, cs_glob_n_ranks);
#else
  snprintf(file_name, 64,
           "partition%cdomain_number_%d",
           _dir_separator, cs_glob_n_ranks);
#endif
  file_name[63] = '\0'; /* Just in case; processor counts would need to be
                           in the exa-range for this to be necessary. */

  /* Test if file exists */

  if (! cs_file_isreg(file_name)) {
    bft_printf(_(" No \"%s\" file available;\n"), file_name);
    if (_use_sfc == false)
      bft_printf(_("   an unoptimized domain partitioning will be used.\n"));
    else
      bft_printf(_("   domain partitioning will use a space-filling curve.\n"));
    return;
  }

  /* Open file */

#if defined(HAVE_MPI)
  rank_pp_in = cs_io_initialize(file_name,
                                "Domain partitioning, R0",
                                CS_IO_MODE_READ,
                                cs_glob_io_hints,
                                echo,
                                cs_glob_mpi_comm);
#else
  rank_pp_in = cs_io_initialize(file_name,
                                "Domain partitioning, R0",
                                CS_IO_MODE_READ,
                                cs_glob_io_hints,
                                echo);
#endif

  if (echo > 0)
    bft_printf("\n");

  /* Loop on read sections */

  while (rank_pp_in != NULL) {

    /* Receive headers */

    cs_io_read_header(rank_pp_in, &header);

    /* Treatment according to the header name */

    if (strncmp(header.sec_name, "n_cells",
                CS_IO_NAME_LEN) == 0) {

      if (header.n_vals != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name,
                  cs_io_get_name(rank_pp_in));
      else {
        cs_io_set_cs_gnum(&header, rank_pp_in);
        cs_io_read_global(&header, &n_g_cells, rank_pp_in);
        if (n_g_cells != mesh->n_g_cells)
          bft_error(__FILE__, __LINE__, 0,
                    _("The number of cells reported by file\n"
                      "\"%s\" (%llu)\n"
                      "does not correspond the those of the mesh (%llu)."),
                    cs_io_get_name(rank_pp_in),
                    (unsigned long long)(n_g_cells),
                    (unsigned long long)(mesh->n_g_cells));
      }

    }
    else if (strncmp(header.sec_name, "n_ranks",
                     CS_IO_NAME_LEN) == 0) {

      if (header.n_vals != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name,
                  cs_io_get_name(rank_pp_in));
      else {
        cs_io_set_cs_lnum(&header, rank_pp_in);
        cs_io_read_global(&header, &n_ranks, rank_pp_in);
        if (n_ranks != cs_glob_n_ranks)
          bft_error(__FILE__, __LINE__, 0,
                    _("The number of ranks reported by file\n"
                      "\"%s\" (%d) does not\n"
                      "correspond to the current number of ranks (%d)."),
                    cs_io_get_name(rank_pp_in), (int)n_ranks,
                    (int)cs_glob_n_ranks);
      }

    }
    else if (strncmp(header.sec_name, "cell:domain number",
                     CS_IO_NAME_LEN) == 0) {

      n_elts = mesh->n_g_cells;
      if (header.n_vals != (cs_file_off_t)n_elts)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name,
                  cs_io_get_name(rank_pp_in));
      else {
        mr->read_cell_rank = 1;
        cs_io_set_cs_lnum(&header, rank_pp_in);
        if (mr->cell_bi.gnum_range[0] > 0)
          n_elts = mr->cell_bi.gnum_range[1] - mr->cell_bi.gnum_range[0];
        BFT_MALLOC(mr->cell_rank, n_elts, cs_lnum_t);
        cs_io_read_block(&header,
                         mr->cell_bi.gnum_range[0],
                         mr->cell_bi.gnum_range[1],
                         mr->cell_rank, rank_pp_in);
        for (i = 0; i < n_elts; i++) /* Convert 1 to n to 0 to n-1 */
          mr->cell_rank[i] -= 1;
      }
      cs_io_finalize(&rank_pp_in);
      rank_pp_in = NULL;
    }

    else
      bft_error(__FILE__, __LINE__, 0,
                _("Section of type <%s> on <%s> is unexpected."),
                header.sec_name, cs_io_get_name(rank_pp_in));
  }

  if (rank_pp_in != NULL)
    cs_io_finalize(&rank_pp_in);
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Mark faces by type (0 for interior, 1 for exterior faces with outwards
 * pointing normal, 2 for exterior faces with inwards pointing normal,
 * 3 for isolated faces) in parallel mode.
 *
 * The mesh structure is also updated with face counts and connectivity sizes.
 *
 * parameters:
 *   mesh              <-> pointer to mesh structure
 *   n_faces           <-- number of local faces
 *   face_ifs          <-- parallel and periodic faces interfaces set
 *   face_cell         <-- local face -> cell connectivity
 *   face_vertices_idx <-- local face -> vertices index
 *   face_type         --> face type marker
 *----------------------------------------------------------------------------*/

static void
_face_type_g(cs_mesh_t                 *mesh,
             cs_lnum_t                  n_faces,
             const cs_interface_set_t  *face_ifs,
             const cs_lnum_t            face_cell[],
             const cs_lnum_t            face_vertices_idx[],
             char                       face_type[])
{
  cs_lnum_t i;
  int j;

  const int n_interfaces = cs_interface_set_size(face_ifs);

  /* Mark base interior faces */

  for (i = 0; i < n_faces; i++) {
    if (face_cell[i*2] > 0 && face_cell[i*2+1] > 0)
      face_type[i] = '\0';
    else if (face_cell[i*2] > 0)
      face_type[i] = '\1';
    else if (face_cell[i*2 + 1] > 0)
      face_type[i] = '\2';
    else {
      face_type[i] = '\3';
    }
  }

  /* Also mark parallel and periodic faces as interior */

  for (j = 0; j < n_interfaces; j++) {

    const cs_interface_t *face_if = cs_interface_set_get(face_ifs, j);
    cs_lnum_t face_if_size = cs_interface_size(face_if);
    const cs_lnum_t *loc_id = cs_interface_get_elt_ids(face_if);

    for (i = 0; i < face_if_size; i++)
      face_type[loc_id[i]] = '\0';

  }

  /* Now count faces of each type */

  mesh->n_i_faces = 0;
  mesh->n_b_faces = 0;
  mesh->i_face_vtx_connect_size = 0;
  mesh->b_face_vtx_connect_size = 0;

  for (i = 0; i < n_faces; i++) {
    cs_lnum_t n_f_vertices = face_vertices_idx[i+1] - face_vertices_idx[i];
    if (face_type[i] == '\0') {
      mesh->n_i_faces += 1;
      mesh->i_face_vtx_connect_size += n_f_vertices;
    }
    else {
      mesh->n_b_faces += 1;
      mesh->b_face_vtx_connect_size += n_f_vertices;
    }
  }
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Mark faces by type (0 for interior, 1 for exterior faces with outwards
 * pointing normal, 2 for exterior faces with inwards pointing normal,
 * 3 for isolated faces) in serial mode.
 *
 * The mesh structure is also updated with face counts and connectivity sizes.
 *
 * parameters:
 *   mesh               <-> pointer to mesh structure
 *   n_faces            <-- number of local faces
 *   n_periodic_couples <-- number of periodic couples associated with
 *                          each periodic list
 *   periodic_couples   <-- array indicating periodic couples (using
 *                          global numberings) for each list
 *   face_cell          <-- local face -> cell connectivity
 *   face_vertices_idx  <-- local face -> vertices index
 *   face_type          --> face type marker
 *----------------------------------------------------------------------------*/

static void
_face_type_l(cs_mesh_t                  *mesh,
             cs_lnum_t                   n_faces,
             const cs_lnum_t             n_periodic_couples[],
             const cs_gnum_t      *const periodic_couples[],
             const cs_lnum_t             face_cell[],
             const cs_lnum_t             face_vertices_idx[],
             char                        face_type[])
{
  cs_lnum_t i;
  int j;

  /* Mark base interior faces */

  for (i = 0; i < n_faces; i++) {
    if (face_cell[i*2] > 0 && face_cell[i*2+1] > 0)
      face_type[i] = '\0';
    else if (face_cell[i*2] > 0)
      face_type[i] = '\1';
    else if (face_cell[i*2 + 1] > 0)
      face_type[i] = '\2';
    else
      face_type[i] = '\3';
  }

  /* Also mark parallel and periodic faces as interior */

  for (i = 0; i < mesh->n_init_perio; i++) {

    const cs_gnum_t *p_couples = periodic_couples[i];

    for (j = 0; j < n_periodic_couples[i]; j++) {
      face_type[p_couples[j*2] - 1] = '\0';
      face_type[p_couples[j*2 + 1] - 1] = '\0';
    }

  }

  /* Now count faces of each type */

  mesh->n_i_faces = 0;
  mesh->n_b_faces = 0;
  mesh->i_face_vtx_connect_size = 0;
  mesh->b_face_vtx_connect_size = 0;

  for (i = 0; i < n_faces; i++) {
    cs_lnum_t n_f_vertices = face_vertices_idx[i+1] - face_vertices_idx[i];
    if (face_type[i] == '\0') {
      mesh->n_i_faces += 1;
      mesh->i_face_vtx_connect_size += n_f_vertices;
    }
    else {
      mesh->n_b_faces += 1;
      mesh->b_face_vtx_connect_size += n_f_vertices;
    }
  }

  mesh->n_g_i_faces = mesh->n_i_faces;
  mesh->n_g_b_faces = mesh->n_b_faces;
}

/*----------------------------------------------------------------------------
 * Build internal and boundary face -> cell connectivity using a common
 * face -> cell connectivity and a face type marker.
 *
 * At this stage, isolated faces, if present, are considered to be
 * boundary faces, as they may participate in future mesh joining
 * operations. Their matching cell number will be set to -1.
 * Remaining isolated faces should be removed before completing
 * the mesh structure.
 *
 * The corresponding arrays in the mesh structure are allocated and
 * defined by this function, and should have been previously empty.
 *
 * parameters:
 *   mesh      <-> pointer to mesh structure
 *   n_faces   <-- number of local faces
 *   face_cell <-- local face -> cell connectivity
 *   face_type <-- face type marker
 *----------------------------------------------------------------------------*/

static void
_extract_face_cell(cs_mesh_t         *mesh,
                   cs_lnum_t          n_faces,
                   const cs_lnum_t    face_cell[],
                   const char         face_type[])
{
  cs_lnum_t i;

  size_t n_i_faces = 0;
  size_t n_b_faces = 0;

  /* Allocate arrays */

  BFT_MALLOC(mesh->i_face_cells, mesh->n_i_faces * 2, cs_int_t);
  BFT_MALLOC(mesh->b_face_cells, mesh->n_b_faces, cs_int_t);

  /* Now copy face -> cell connectivity */

  for (i = 0; i < n_faces; i++) {

    if (face_type[i] == '\0') {
      mesh->i_face_cells[n_i_faces*2]     = face_cell[i*2];
      mesh->i_face_cells[n_i_faces*2 + 1] = face_cell[i*2 + 1];
      n_i_faces++;
    }

    else if (face_type[i] == '\1') {
      mesh->b_face_cells[n_b_faces] = face_cell[i*2];
      n_b_faces++;
    }

    else if (face_type[i] == '\2') {
      mesh->b_face_cells[n_b_faces] = face_cell[i*2 + 1];
      n_b_faces++;
    }

    else if (face_type[i] == '\3') {
      mesh->b_face_cells[n_b_faces] = -1;
      mesh->n_g_free_faces += 1;
      n_b_faces++;
    }
  }
}

/*----------------------------------------------------------------------------
 * Build internal and boundary face -> vertices connectivity using a common
 * face -> vertices connectivity and a face type marker.
 *
 * The corresponding arrays in the mesh structure are allocated and
 * defined by this function, and should have been previously empty.
 *
 * parameters:
 *   mesh              <-> pointer to mesh structure
 *   n_faces           <-- number of local faces
 *   face_vertices_idx <-- local face -> vertices index
 *   face_vertices     <-- local face -> vertices connectivity
 *   face_type         <-- face type marker
 *----------------------------------------------------------------------------*/

static void
_extract_face_vertices(cs_mesh_t         *mesh,
                       cs_lnum_t          n_faces,
                       const cs_lnum_t    face_vertices_idx[],
                       const cs_lnum_t    face_vertices[],
                       const char         face_type[])
{
  cs_lnum_t i;
  size_t j;

  size_t n_i_faces = 0;
  size_t n_b_faces = 0;

  /* Allocate and initialize */

  BFT_MALLOC(mesh->i_face_vtx_idx, mesh->n_i_faces+1, cs_int_t);
  BFT_MALLOC(mesh->i_face_vtx_lst, mesh->i_face_vtx_connect_size, cs_int_t);

  mesh->i_face_vtx_idx[0] = 1;

  BFT_MALLOC(mesh->b_face_vtx_idx, mesh->n_b_faces+1, cs_int_t);
  mesh->b_face_vtx_idx[0] = 1;

  if (mesh->n_b_faces > 0)
    BFT_MALLOC(mesh->b_face_vtx_lst, mesh->b_face_vtx_connect_size, cs_int_t);

  /* Now copy face -> vertices connectivity */

  for (i = 0; i < n_faces; i++) {

    size_t n_f_vertices = face_vertices_idx[i+1] - face_vertices_idx[i];
    const cs_lnum_t *_face_vtx = face_vertices + face_vertices_idx[i];

    if (face_type[i] == '\0') {
      cs_lnum_t *_i_face_vtx =   mesh->i_face_vtx_lst
                                + mesh->i_face_vtx_idx[n_i_faces] - 1;
      for (j = 0; j < n_f_vertices; j++)
        _i_face_vtx[j] = _face_vtx[j];
      mesh->i_face_vtx_idx[n_i_faces + 1] =   mesh->i_face_vtx_idx[n_i_faces]
                                            + n_f_vertices;
      n_i_faces++;
    }

    else if (face_type[i] == '\1' || face_type[i] == '\3') {
      cs_lnum_t *_b_face_vtx =   mesh->b_face_vtx_lst
                                + mesh->b_face_vtx_idx[n_b_faces] - 1;
      for (j = 0; j < n_f_vertices; j++)
        _b_face_vtx[j] = _face_vtx[j];
      mesh->b_face_vtx_idx[n_b_faces + 1] =   mesh->b_face_vtx_idx[n_b_faces]
                                            + n_f_vertices;
      n_b_faces++;
    }

    else if (face_type[i] == '\2') {
      cs_lnum_t *_b_face_vtx =   mesh->b_face_vtx_lst
                                + mesh->b_face_vtx_idx[n_b_faces] - 1;
      for (j = 0; j < n_f_vertices; j++)
        _b_face_vtx[j] = _face_vtx[n_f_vertices - j - 1];
      mesh->b_face_vtx_idx[n_b_faces + 1] =   mesh->b_face_vtx_idx[n_b_faces]
                                            + n_f_vertices;
      n_b_faces++;
    }

  }
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Build internal and boundary face -> global numberings using a common
 * face group class id and a face type marker.
 *
 * The corresponding arrays in the mesh structure are allocated and
 * defined by this function, and should have been previously empty.
 *
 * parameters:
 *   mesh            <-> pointer to mesh structure
 *   n_faces         <-- number of local faces
 *   global_face_num <-- global face numbers
 *   face_type       <-- face type marker
 *----------------------------------------------------------------------------*/

static void
_extract_face_gnum(cs_mesh_t         *mesh,
                   cs_lnum_t          n_faces,
                   const cs_gnum_t    global_face_num[],
                   const char         face_type[])
{
  cs_lnum_t i;

  size_t n_i_faces = 0;
  size_t n_b_faces = 0;

  cs_lnum_t *global_i_face = NULL;
  cs_lnum_t *global_b_face = NULL;

  fvm_io_num_t *tmp_face_num = NULL;

  /* Allocate arrays (including temporary arrays) */

  BFT_MALLOC(mesh->global_i_face_num, mesh->n_i_faces, cs_gnum_t);
  BFT_MALLOC(mesh->global_b_face_num, mesh->n_b_faces, cs_gnum_t);

  BFT_MALLOC(global_i_face, mesh->n_i_faces, cs_lnum_t);
  BFT_MALLOC(global_b_face, mesh->n_b_faces, cs_lnum_t);

  /* Now build internal and boundary face lists */

  for (i = 0; i < n_faces; i++) {

    if (face_type[i] == '\0')
      global_i_face[n_i_faces++] = i+1;

    else
      global_b_face[n_b_faces++] = i+1;

  }

  /* Build an I/O numbering on internal faces to compact the global numbering */

  tmp_face_num = fvm_io_num_create(global_i_face,
                                   global_face_num,
                                   n_i_faces,
                                   0);

  memcpy(mesh->global_i_face_num,
         fvm_io_num_get_global_num(tmp_face_num),
         n_i_faces*sizeof(cs_gnum_t));

  mesh->n_g_i_faces = fvm_io_num_get_global_count(tmp_face_num);

  assert(fvm_io_num_get_local_count(tmp_face_num) == (cs_lnum_t)n_i_faces);

  tmp_face_num = fvm_io_num_destroy(tmp_face_num);

  /* Build an I/O numbering on boundary faces to compact the global numbering */

  tmp_face_num = fvm_io_num_create(global_b_face,
                                   global_face_num,
                                   n_b_faces,
                                   0);

  if (n_b_faces > 0)
    memcpy(mesh->global_b_face_num,
           fvm_io_num_get_global_num(tmp_face_num),
           n_b_faces*sizeof(cs_gnum_t));

  mesh->n_g_b_faces = fvm_io_num_get_global_count(tmp_face_num);

  assert(fvm_io_num_get_local_count(tmp_face_num) == (cs_lnum_t)n_b_faces);

  tmp_face_num = fvm_io_num_destroy(tmp_face_num);

  /* Free remaining temporary arrays */

  BFT_FREE(global_i_face);
  BFT_FREE(global_b_face);
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Build internal and boundary face -> group class id using a common
 * face group class id and a face type marker.
 *
 * The corresponding arrays in the mesh structure are allocated and
 * defined by this function, and should have been previously empty.
 *
 * parameters:
 *   mesh       <-> pointer to mesh structure
 *   n_faces    <-- number of local faces
 *   face_gc_id <-- local face group class id
 *   face_type  <-- face type marker
 *----------------------------------------------------------------------------*/

static void
_extract_face_gc_id(cs_mesh_t        *mesh,
                   cs_lnum_t          n_faces,
                   const cs_lnum_t    face_gc_id[],
                   const char         face_type[])
{
  cs_lnum_t i;

  size_t n_i_faces = 0;
  size_t n_b_faces = 0;

  /* Allocate arrays */

  BFT_MALLOC(mesh->i_face_family, mesh->n_i_faces, cs_int_t);
  BFT_MALLOC(mesh->b_face_family, mesh->n_b_faces, cs_int_t);

  /* Now copy face group class (family) id */

  for (i = 0; i < n_faces; i++) {

    assert(face_gc_id[i] > -1 && face_gc_id[i] <= mesh->n_families);

    if (face_type[i] == '\0')
      mesh->i_face_family[n_i_faces++] = face_gc_id[i];

    else
      mesh->b_face_family[n_b_faces++] = face_gc_id[i];

  }
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Renumber face interface references from mixed faces to interior faces.
 *
 * parameters:
 *   face_ifs          <-> parallel and periodic faces interfaces set
 *   n_faces           <-- number of local faces
 *   face_type         <-- face type marker
 *----------------------------------------------------------------------------*/

static void
_face_ifs_to_interior(cs_interface_set_t  *face_ifs,
                      cs_lnum_t            n_faces,
                      const char           face_type[])
{
  cs_lnum_t i;

  cs_lnum_t   i_face_count = 0;
  cs_lnum_t  *i_face_id = NULL;

  /* Build face renumbering */

  BFT_MALLOC(i_face_id, n_faces, cs_lnum_t);

  for (i = 0; i < n_faces; i++) {
    if (face_type[i] == '\0')
      i_face_id[i] = i_face_count++;
    else
      i_face_id[i] = -1;
  }

  cs_interface_set_renumber(face_ifs, i_face_id);

  BFT_FREE(i_face_id);
}

/*----------------------------------------------------------------------------
 * Compare periodic couples in global numbering form (qsort function).
 *
 * parameters:
 *   x <-> pointer to first couple
 *   y <-> pointer to second couple
 *
 * returns:
 *   lexicographical
 *----------------------------------------------------------------------------*/

static int _compare_couples(const void *x, const void *y)
{
  int retval = 1;

  const cs_gnum_t *c0 = x;
  const cs_gnum_t *c1 = y;

  if (c0[0] < c1[0])
    retval = -1;

  else if (c0[0] == c1[0]) {
    if (c0[1] < c1[1])
      retval = -1;
    else if (c0[1] == c1[1])
      retval = 0;
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Extract periodic face connectivity information for mesh builder when
 * running in parallel mode.
 *
 * parameters:
 *   mesh     <-- pointer to mesh structure
 *   mb       <-> pointer to mesh builder structure
 *   face_ifs <-- parallel and periodic faces interfaces set
 *----------------------------------------------------------------------------*/

static void
_extract_periodic_faces_g(const cs_mesh_t           *mesh,
                          cs_mesh_builder_t         *mb,
                          const cs_interface_set_t  *face_ifs)
{
  int i, j;
  cs_lnum_t k, l;

  int perio_count = 0;
  cs_lnum_t  *send_index = NULL;
  cs_gnum_t  *recv_num = NULL;
  int  *tr_id = NULL;

  cs_datatype_t gnum_type = CS_GNUM_TYPE;

  const int n_perio = mesh->n_init_perio;
  const int n_interfaces = cs_interface_set_size(face_ifs);
  const cs_gnum_t *face_gnum = mesh->global_i_face_num;

  /* Allocate arrays in mesh builder (initializing per_face_idx) */

  assert(mesh->periodicity != NULL);
  assert(mb != NULL);
  assert(mb->n_perio == 0);

  mb->n_perio = n_perio;

  BFT_MALLOC(mb->n_perio_couples, n_perio, cs_lnum_t);
  BFT_MALLOC(mb->perio_couples, n_perio, cs_gnum_t *);

  for (i = 0; i < n_perio; i++) {
    mb->n_perio_couples[i] = 0;
    mb->perio_couples[i] = NULL;
  }

  /* List direct and reverse transforms */

  BFT_MALLOC(tr_id, n_perio*2, int);

  for (i = 0; i < n_perio*2; i++) {
    int rev_id = fvm_periodicity_get_reverse_id(mesh->periodicity, i);
    if (i < rev_id) {
      int parent_ids[2];
      fvm_periodicity_get_parent_ids(mesh->periodicity, i, parent_ids);
      if (parent_ids[0] < 0 && parent_ids[1] < 0) {
        tr_id[perio_count*2] = i + 1;
        tr_id[perio_count*2 + 1] = rev_id + 1;
        perio_count++;
      }
    }
  }
  assert(perio_count == n_perio);

  for (i = 0; i < n_interfaces; i++) {
    const cs_interface_t *face_if = cs_interface_set_get(face_ifs, i);
    const cs_lnum_t *tr_index = cs_interface_get_tr_index(face_if);
    for (j = 0; j < n_perio; j++) {
      const cs_lnum_t n_tr_faces = (  tr_index[tr_id[j*2] + 1]
                                    - tr_index[tr_id[j*2]]);
      mb->n_perio_couples[j] += n_tr_faces;
    }
  }

  BFT_MALLOC(recv_num, cs_interface_set_n_elts(face_ifs), cs_gnum_t);

  cs_interface_set_copy_array(face_ifs,
                              gnum_type,
                              1,
                              true, /* src_on_parent */
                              face_gnum,
                              recv_num);

  /* Prepare send buffer (send reverse transformation values) */

  BFT_FREE(send_index);

  for (i = 0; i < n_perio; i++)
    BFT_MALLOC(mb->perio_couples[i], mb->n_perio_couples[i]*2, cs_gnum_t);

  /* Reset couples count */

  for (i = 0; i < n_perio; i++)
    mb->n_perio_couples[i] = 0;

  /* Copy face couples to mesh builder */

  for (i = 0, j = 0, l = 0; i < n_interfaces; i++) {

    const cs_interface_t *face_if = cs_interface_set_get(face_ifs, i);
    const cs_lnum_t *tr_index = cs_interface_get_tr_index(face_if);
    const cs_lnum_t *elt_id = cs_interface_get_elt_ids(face_if);

    l += tr_index[1];

    for (j = 0; j < n_perio; j++) {

      /* Count couples in direct periodicity */

      cs_lnum_t nc = mb->n_perio_couples[j]*2;
      const cs_lnum_t start_id = tr_index[tr_id[j*2]];
      const cs_lnum_t end_id = tr_index[tr_id[j*2] + 1];

      for (k = start_id; k < end_id; k++) {
        cs_lnum_t f_id = elt_id[k];
        mb->perio_couples[j][nc++] = face_gnum[f_id];
        mb->perio_couples[j][nc++] = recv_num[l++];
      }
      mb->n_perio_couples[j] = nc/2;

      /* Ignore couples in reverse periodicity */

      l += tr_index[tr_id[j*2 + 1] + 1] - tr_index[tr_id[j*2 + 1]];

    }

  }

  BFT_FREE(recv_num);
  BFT_FREE(tr_id);

  /* Now sort couples in place for future use (more for consistency
     and ease of verification than absolutely necessary) */

  for (i = 0; i < n_perio; i++) {
    if (mb->n_perio_couples[i] > 0)
      qsort(mb->perio_couples[i],
            mb->n_perio_couples[i],
            sizeof(cs_gnum_t) * 2,
            &_compare_couples);
  }
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Extract periodic face connectivity information for mesh builder when
 * running in serial mode.
 *
 * Arrays are simply transferred from the mesh reader to the builder and
 * renumbered
 *
 * parameters:
 *   mr            <-> pointer to mesh reader structure
 *   mb            <-> pointer to mesh builder structure
 *   n_init_perio  <-- number of initial periodicities
 *   n_faces       <-- number of local faces
 *   face_type     <-- face type marker
 *----------------------------------------------------------------------------*/

static void
_extract_periodic_faces_l(_mesh_reader_t     *mr,
                          cs_mesh_builder_t  *mb,
                          int                 n_init_perio,
                          cs_lnum_t           n_faces,
                          const char          face_type[])
{
  int i;

  cs_gnum_t   next_face_num = 1;
  cs_gnum_t  *i_face_num = NULL;

  /* Transfer arrays from reader to builder, then renumber couples */

  assert(mb != NULL);
  assert(mb->n_perio_couples == NULL);

  mb->n_perio = n_init_perio;
  mb->n_perio_couples = mr->n_per_face_couples;
  mb->perio_couples = mr->per_face_couples;

  mr->n_per_face_couples = NULL;
  mr->per_face_couples = NULL;

  /* Build face renumbering */

  BFT_MALLOC(i_face_num, n_faces, cs_gnum_t);

  for (i = 0; i < n_faces; i++) {
    if (face_type[i] == '\0')
      i_face_num[i] = next_face_num++;
    else
      i_face_num[i] = 0;
  }

  /* Apply new numbering */

  for (i = 0; i < n_init_perio; i++) {

    size_t j;
    cs_gnum_t *p_couples = mb->perio_couples[i];
    const size_t n_vals = mb->n_perio_couples[i] * 2;

    for (j = 0; j < n_vals; j++) {
      p_couples[j] = i_face_num[p_couples[j] - 1];
      assert(p_couples[j] > 0);
    }
  }

  BFT_FREE(i_face_num);
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Compute cell centers using minimal local data.
 *
 * parameters:
 *   n_cells      <-- number of cells
 *   n_faces      <-- number of faces
 *   face_cells   <-- face -> cells connectivity
 *   face_vtx_idx <-- face -> vertices connectivity index
 *   face_vtx     <-- face -> vertices connectivity
 *   vtx_coord    <-- vertex coordinates
 *   cell_center  --> cell centers
 *----------------------------------------------------------------------------*/

static void
_cell_center(cs_lnum_t         n_cells,
             cs_lnum_t         n_faces,
             const cs_lnum_t   face_cells[],
             const cs_lnum_t   face_vtx_idx[],
             const cs_lnum_t   face_vtx[],
             const cs_real_t   vtx_coord[],
             cs_coord_t        cell_center[])
{
  cs_lnum_t i, j;
  cs_lnum_t vtx_id, face_id, start_id, end_id;
  cs_lnum_t n_face_vertices;
  cs_coord_t ref_normal[3], vtx_cog[3], face_center[3];

  cs_lnum_t n_max_face_vertices = 0;

  _vtx_coords_t *face_vtx_coord = NULL;
  cs_coord_t *weight = NULL;

  const double surf_epsilon = 1e-24;

  assert(face_vtx_idx[0] == 0);

  BFT_MALLOC(weight, n_cells, cs_coord_t);

  for (i = 0; i < n_cells; i++) {
    weight[i] = 0.0;
    for (j = 0; j < 3; j++)
      cell_center[i*3 + j] = 0.0;
  }

  /* Counting and allocation */

  n_max_face_vertices = 0;

  for (face_id = 0; face_id < n_faces; face_id++) {
    n_face_vertices = face_vtx_idx[face_id + 1] - face_vtx_idx[face_id];
    if (n_max_face_vertices <= n_face_vertices)
      n_max_face_vertices = n_face_vertices;
  }

  BFT_MALLOC(face_vtx_coord, n_max_face_vertices, _vtx_coords_t);

  /* Loop on each face */

  for (face_id = 0; face_id < n_faces; face_id++) {

    /* Initialization */

    cs_lnum_t tri_id;

    cs_lnum_t cell_id_0 = face_cells[face_id*2] -1;
    cs_lnum_t cell_id_1 = face_cells[face_id*2 + 1] -1;
    cs_coord_t unweighted_center[3] = {0.0, 0.0, 0.0};
    cs_coord_t face_surface = 0.0;

    n_face_vertices = 0;

    start_id = face_vtx_idx[face_id];
    end_id = face_vtx_idx[face_id + 1];

    /* Define the polygon (P) according to the vertices (Pi) of the face */

    for (vtx_id = start_id; vtx_id < end_id; vtx_id++) {

      cs_lnum_t shift = 3 * (face_vtx[vtx_id] - 1);
      for (i = 0; i < 3; i++)
        face_vtx_coord[n_face_vertices][i] = vtx_coord[shift + i];
      n_face_vertices++;

    }

    /* Compute the barycentre of the face vertices */

    for (i = 0; i < 3; i++) {
      vtx_cog[i] = 0.0;
      for (vtx_id = 0; vtx_id < n_face_vertices; vtx_id++)
        vtx_cog[i] += face_vtx_coord[vtx_id][i];
      vtx_cog[i] /= n_face_vertices;
    }

    /* Loop on the triangles of the face (defined by an edge of the face
       and its barycentre) */

    for (i = 0; i < 3; i++) {
      ref_normal[i] = 0.;
      face_center[i] = 0.0;
    }

    for (tri_id = 0 ; tri_id < n_face_vertices ; tri_id++) {

      cs_coord_t tri_surface;
      cs_coord_t vect1[3], vect2[3], tri_normal[3], tri_center[3];

      cs_lnum_t id0 = tri_id;
      cs_lnum_t id1 = (tri_id + 1)%n_face_vertices;

      /* Normal for each triangle */

      for (i = 0; i < 3; i++) {
        vect1[i] = face_vtx_coord[id0][i] - vtx_cog[i];
        vect2[i] = face_vtx_coord[id1][i] - vtx_cog[i];
      }

      tri_normal[0] = vect1[1] * vect2[2] - vect2[1] * vect1[2];
      tri_normal[1] = vect2[0] * vect1[2] - vect1[0] * vect2[2];
      tri_normal[2] = vect1[0] * vect2[1] - vect2[0] * vect1[1];

      if (tri_id == 0) {
        for (i = 0; i < 3; i++)
          ref_normal[i] = tri_normal[i];
      }

      /* Center of gravity for a triangle */

      for (i = 0; i < 3; i++) {
        tri_center[i] = (  vtx_cog[i]
                         + face_vtx_coord[id0][i]
                         + face_vtx_coord[id1][i]) / 3.0;
      }

      tri_surface = sqrt(  tri_normal[0]*tri_normal[0]
                         + tri_normal[1]*tri_normal[1]
                         + tri_normal[2]*tri_normal[2]) * 0.5;

      if ((  tri_normal[0]*ref_normal[0]
           + tri_normal[1]*ref_normal[1]
           + tri_normal[2]*ref_normal[2]) < 0.0)
        tri_surface *= -1.0;

      /* Now compute contribution to face center and surface */

      face_surface += tri_surface;

      for (i = 0; i < 3; i++) {
        face_center[i] += tri_surface * tri_center[i];
        unweighted_center[i] = tri_center[i];
      }

    } /* End of loop  on triangles of the face */

    if (face_surface > surf_epsilon) {
      for (i = 0; i < 3; i++)
        face_center[i] /= face_surface;
    }
    else {
      face_surface = surf_epsilon;
for (i = 0; i < 3; i++)
        face_center[i] = unweighted_center[i] * face_surface / n_face_vertices;
    }

    /* Now contribute to cell centers */

    assert(cell_id_0 > -2 && cell_id_1 > -2);

    if (cell_id_0 > -1) {
      for (i = 0; i < 3; i++)
        cell_center[cell_id_0*3 + i] += face_center[i]*face_surface;
      weight[cell_id_0] += face_surface;
    }

    if (cell_id_1 > -1) {
      for (i = 0; i < 3; i++)
        cell_center[cell_id_1*3 + i] += face_center[i]*face_surface;
      weight[cell_id_1] += face_surface;
    }

  } /* End of loop on faces */

  BFT_FREE(face_vtx_coord);

  for (i = 0; i < n_cells; i++) {
    for (j = 0; j < 3; j++)
      cell_center[i*3 + j] /= weight[i];
  }

  BFT_FREE(weight);
}

/*----------------------------------------------------------------------------
 * Compute cell centers using block data read from file.
 *
 * parameters:
 *   mr          <-- pointer to mesh reader helper structure
 *   cell_center --> cell centers array
 *   comm        <-- associated MPI communicator
 *----------------------------------------------------------------------------*/

static void
_precompute_cell_center(const _mesh_reader_t    *mr,
                        cs_coord_t               cell_center[],
                        MPI_Comm                 comm)
{
  cs_lnum_t i;
  int n_ranks = 0;

  cs_datatype_t gnum_type = (sizeof(cs_gnum_t) == 8) ? CS_UINT64 : CS_UINT32;
  cs_datatype_t real_type = (sizeof(cs_real_t) == 8) ? CS_DOUBLE : CS_FLOAT;

  cs_lnum_t _n_cells = 0;
  cs_lnum_t _n_faces = 0;
  cs_lnum_t _n_vertices = 0;

  cs_gnum_t *_cell_num = NULL;
  cs_gnum_t *_face_num = NULL;
  cs_gnum_t *_vtx_num = NULL;
  cs_gnum_t *_face_gcells = NULL;
  cs_gnum_t *_face_gvertices = NULL;

  cs_lnum_t *_face_cells = NULL;
  cs_lnum_t *_face_vertices_idx = NULL;
  cs_lnum_t *_face_vertices = NULL;

  cs_real_t *_vtx_coord = NULL;

  cs_block_to_part_t *d = NULL;

  /* Initialization */

  MPI_Comm_size(comm, &n_ranks);

  assert((sizeof(cs_lnum_t) == 4) || (sizeof(cs_lnum_t) == 8));

  _n_cells = mr->cell_bi.gnum_range[1] - mr->cell_bi.gnum_range[0];

  BFT_MALLOC(_cell_num, _n_cells, cs_gnum_t);

  for (i = 0; i < _n_cells; i++)
    _cell_num[i] = mr->cell_bi.gnum_range[0] + i;

  if (_n_cells == 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Number of cells on rank %d is zero.\n"
                "(number of cells / number of processes ratio too low)."),
              (int)cs_glob_rank_id);

  /* Distribute faces */
  /*------------------*/

  d = cs_block_to_part_create_by_adj_s(comm,
                                       mr->face_bi,
                                       mr->cell_bi,
                                       2,
                                       mr->face_cells,
                                       NULL,
                                       NULL);

  _n_faces = cs_block_to_part_get_n_part_ents(d);

  BFT_MALLOC(_face_gcells, _n_faces*2, cs_gnum_t);

  /* Face -> cell connectivity */

  cs_block_to_part_copy_array(d,
                              gnum_type,
                              2,
                              mr->face_cells,
                              _face_gcells);

  /* Now convert face -> cell connectivity to local cell numbers */

  BFT_MALLOC(_face_cells, _n_faces*2, cs_lnum_t);

  cs_block_to_part_global_to_local(_n_faces*2,
                                   1,
                                   _n_cells,
                                   _cell_num,
                                   _face_gcells,
                                   _face_cells);

  BFT_FREE(_cell_num);
  BFT_FREE(_face_gcells);

  /* Face connectivity */

  BFT_MALLOC(_face_vertices_idx, _n_faces + 1, cs_lnum_t);

  cs_block_to_part_copy_index(d,
                              mr->face_vertices_idx,
                              _face_vertices_idx);

  BFT_MALLOC(_face_gvertices, _face_vertices_idx[_n_faces], cs_gnum_t);

  cs_block_to_part_copy_indexed(d,
                                gnum_type,
                                mr->face_vertices_idx,
                                mr->face_vertices,
                                _face_vertices_idx,
                                _face_gvertices);

  _face_num = cs_block_to_part_transfer_gnum(d);

  cs_block_to_part_destroy(&d);

  /* Vertices */

  d = cs_block_to_part_create_adj(comm,
                                  mr->vertex_bi,
                                  _face_vertices_idx[_n_faces],
                                  _face_gvertices);

  _n_vertices = cs_block_to_part_get_n_part_ents(d);

  BFT_MALLOC(_vtx_coord, _n_vertices*3, cs_real_t);

  cs_block_to_part_copy_array(d,
                              real_type,
                              3,
                              mr->vertex_coords,
                              _vtx_coord);

  _vtx_num = cs_block_to_part_transfer_gnum(d);

  cs_block_to_part_destroy(&d);

  /* Now convert face -> vertex connectivity to local vertex numbers */

  BFT_MALLOC(_face_vertices, _face_vertices_idx[_n_faces], cs_lnum_t);

  cs_block_to_part_global_to_local(_face_vertices_idx[_n_faces],
                                   1,
                                   _n_vertices,
                                   _vtx_num,
                                   _face_gvertices,
                                   _face_vertices);

  BFT_FREE(_face_gvertices);

  _cell_center(_n_cells,
               _n_faces,
               _face_cells,
               _face_vertices_idx,
               _face_vertices,
               _vtx_coord,
               cell_center);

  BFT_FREE(_vtx_coord);
  BFT_FREE(_vtx_num);

  BFT_FREE(_face_cells);

  BFT_FREE(_face_vertices_idx);
  BFT_FREE(_face_vertices);

  BFT_FREE(_face_num);
}

/*----------------------------------------------------------------------------
 * Compute free (isolated) face centers using minimal local data.
 *
 * parameters:
 *   n_f_faces     <-- number of faces
 *   f_face_ids    <-- list of free faces
 *   face_vtx_idx  <-- face -> vertices connectivity index
 *   face_vtx      <-- face -> vertices connectivity
 *   vtx_coord     <-- vertex coordinates
 *   f_face_center --> free face centers
 *----------------------------------------------------------------------------*/

static void
_f_face_center(cs_lnum_t         n_f_faces,
               cs_lnum_t         f_face_ids[],
               const cs_lnum_t   face_vtx_idx[],
               const cs_lnum_t   face_vtx[],
               const cs_real_t   vtx_coord[],
               cs_coord_t        f_face_center[])
{
  cs_lnum_t i, j, k;
  cs_lnum_t vtx_id, start_id, end_id;
  cs_lnum_t n_face_vertices;
  cs_coord_t ref_normal[3], vtx_cog[3];

  cs_lnum_t n_max_face_vertices = 0;

  _vtx_coords_t *face_vtx_coord = NULL;

  const double surf_epsilon = 1e-24;

  assert(face_vtx_idx[0] == 0);

  for (i = 0; i < n_f_faces; i++) {
    for (j = 0; j < 3; j++)
      f_face_center[i*3 + j] = 0.0;
  }

  /* Counting and allocation */

  n_max_face_vertices = 0;

  for (k = 0; k < n_f_faces; k++) {
    cs_lnum_t face_id = f_face_ids[k];
    n_face_vertices = face_vtx_idx[face_id + 1] - face_vtx_idx[face_id];
    if (n_max_face_vertices <= n_face_vertices)
      n_max_face_vertices = n_face_vertices;
  }

  BFT_MALLOC(face_vtx_coord, n_max_face_vertices, _vtx_coords_t);

  /* Loop on each face */

  for (k = 0; k < n_f_faces; k++) {

    cs_lnum_t tri_id;

    /* Initialization */

    cs_lnum_t face_id = f_face_ids[k];
    cs_coord_t unweighted_center[3] = {0.0, 0.0, 0.0};
    cs_coord_t face_surface = 0.0;
    cs_coord_t *face_center = f_face_center + (k*3);

    n_face_vertices = 0;

    start_id = face_vtx_idx[face_id];
    end_id = face_vtx_idx[face_id + 1];

    /* Define the polygon (P) according to the vertices (Pi) of the face */

    for (vtx_id = start_id; vtx_id < end_id; vtx_id++) {

      cs_lnum_t shift = 3 * (face_vtx[vtx_id] - 1);
      for (i = 0; i < 3; i++)
        face_vtx_coord[n_face_vertices][i] = vtx_coord[shift + i];
      n_face_vertices++;

    }

    /* Compute the barycentre of the face vertices */

    for (i = 0; i < 3; i++) {
      vtx_cog[i] = 0.0;
      for (vtx_id = 0; vtx_id < n_face_vertices; vtx_id++)
        vtx_cog[i] += face_vtx_coord[vtx_id][i];
      vtx_cog[i] /= n_face_vertices;
    }

    /* Loop on the triangles of the face (defined by an edge of the face
       and its barycentre) */

    for (i = 0; i < 3; i++) {
      ref_normal[i] = 0.;
      face_center[i] = 0.0;
    }

    for (tri_id = 0 ; tri_id < n_face_vertices ; tri_id++) {

      cs_coord_t tri_surface;
      cs_coord_t vect1[3], vect2[3], tri_normal[3], tri_center[3];

      cs_lnum_t id0 = tri_id;
      cs_lnum_t id1 = (tri_id + 1)%n_face_vertices;

      /* Normal for each triangle */

      for (i = 0; i < 3; i++) {
        vect1[i] = face_vtx_coord[id0][i] - vtx_cog[i];
        vect2[i] = face_vtx_coord[id1][i] - vtx_cog[i];
      }

      tri_normal[0] = vect1[1] * vect2[2] - vect2[1] * vect1[2];
      tri_normal[1] = vect2[0] * vect1[2] - vect1[0] * vect2[2];
      tri_normal[2] = vect1[0] * vect2[1] - vect2[0] * vect1[1];

      if (tri_id == 0) {
        for (i = 0; i < 3; i++)
          ref_normal[i] = tri_normal[i];
      }

      /* Center of gravity for a triangle */

      for (i = 0; i < 3; i++) {
        tri_center[i] = (  vtx_cog[i]
                         + face_vtx_coord[id0][i]
                         + face_vtx_coord[id1][i]) / 3.0;
      }

      tri_surface = sqrt(  tri_normal[0]*tri_normal[0]
                         + tri_normal[1]*tri_normal[1]
                         + tri_normal[2]*tri_normal[2]) * 0.5;

      if ((  tri_normal[0]*ref_normal[0]
           + tri_normal[1]*ref_normal[1]
           + tri_normal[2]*ref_normal[2]) < 0.0)
        tri_surface *= -1.0;

      /* Now compute contribution to face center and surface */

      face_surface += tri_surface;

      for (i = 0; i < 3; i++) {
        face_center[i] += tri_surface * tri_center[i];
        unweighted_center[i] = tri_center[i];
      }

    } /* End of loop  on triangles of the face */

    if (face_surface > surf_epsilon) {
      for (i = 0; i < 3; i++)
        face_center[i] /= face_surface;
    }
    else {
      face_surface = surf_epsilon;
      for (i = 0; i < 3; i++)
        face_center[i] = unweighted_center[i] * face_surface / n_face_vertices;
    }
  } /* End of loop on faces */

  BFT_FREE(face_vtx_coord);
}

/*----------------------------------------------------------------------------
 * Compute face centers using block data read from file.
 *
 * parameters:
 *   mr          <-- pointer to mesh reader helper structure
 *   n_f_faces   <-- local number of free faces
 *   f_face_ids  <-- free face ids
 *   face_center --> cell centers array
 *   comm        <-- associated MPI communicator
 *----------------------------------------------------------------------------*/

static void
_precompute_free_face_center(const _mesh_reader_t   *mr,
                             cs_lnum_t               n_f_faces,
                             cs_lnum_t               f_face_ids[],
                             cs_coord_t              f_face_center[],
                             MPI_Comm                comm)
{
  int n_ranks = 0;

  cs_datatype_t real_type = (sizeof(cs_real_t) == 8) ? CS_DOUBLE : CS_FLOAT;

  cs_lnum_t _n_faces = 0;
  cs_lnum_t _n_vertices = 0;

  cs_gnum_t *_vtx_num = NULL;
  cs_lnum_t *_face_vertices = NULL;

  cs_real_t *_vtx_coord = NULL;

  cs_block_to_part_t *d = NULL;

  /* Initialization */

  MPI_Comm_size(comm, &n_ranks);

  assert((sizeof(cs_lnum_t) == 4) || (sizeof(cs_lnum_t) == 8));

  _n_faces = mr->face_bi.gnum_range[1] - mr->face_bi.gnum_range[0];

  /* Distribute vertices */
  /*---------------------*/

  d = cs_block_to_part_create_adj(comm,
                                  mr->vertex_bi,
                                  mr->face_vertices_idx[_n_faces],
                                  mr->face_vertices);

  _n_vertices = cs_block_to_part_get_n_part_ents(d);

  BFT_MALLOC(_vtx_coord, _n_vertices*3, cs_real_t);

  cs_block_to_part_copy_array(d,
                               real_type,
                               3,
                               mr->vertex_coords,
                               _vtx_coord);

  _vtx_num = cs_block_to_part_transfer_gnum(d);

  cs_block_to_part_destroy(&d);

  /* Now convert face -> vertex connectivity to local vertex numbers */

  BFT_MALLOC(_face_vertices, mr->face_vertices_idx[_n_faces], cs_lnum_t);

  cs_block_to_part_global_to_local(mr->face_vertices_idx[_n_faces],
                                   1,
                                   _n_vertices,
                                   _vtx_num,
                                   mr->face_vertices,
                                   _face_vertices);

  _f_face_center(n_f_faces,
                 f_face_ids,
                 mr->face_vertices_idx,
                 _face_vertices,
                 _vtx_coord,
                 f_face_center);

  BFT_FREE(_vtx_coord);
  BFT_FREE(_vtx_num);
  BFT_FREE(_face_vertices);
}

/*----------------------------------------------------------------------------
 * Compute cell centers using block data read from file.
 *
 * parameters:
 *   mr          <-_ pointer to mesh reader helper structure
 *   cell_rank   --> cell rank
 *   comm        <-- associated MPI communicator
 *----------------------------------------------------------------------------*/

static void
_cell_rank_by_sfc(const _mesh_reader_t     *mr,
                  int                       cell_rank[],
                  MPI_Comm                  comm)
{
  cs_lnum_t i;
  cs_lnum_t n_cells = 0, block_size = 0, rank_step = 0;
  cs_coord_t *cell_center = NULL;
  fvm_io_num_t *cell_io_num = NULL;
  const cs_gnum_t *cell_num = NULL;

  bft_printf(_(" Partitioning by space-filling curve: %s.\n"),
             fvm_io_num_sfc_type_name[_sfc_type]);

  n_cells = mr->cell_bi.gnum_range[1] - mr->cell_bi.gnum_range[0];
  block_size = mr->cell_bi.block_size;
  rank_step = mr->cell_bi.rank_step;

  BFT_MALLOC(cell_center, n_cells*3, cs_coord_t);

  _precompute_cell_center(mr, cell_center, comm);

  cell_io_num = fvm_io_num_create_from_sfc(cell_center,
                                           3,
                                           n_cells,
                                           _sfc_type);

  BFT_FREE(cell_center);

  cell_num = fvm_io_num_get_global_num(cell_io_num);

  /* Determine rank based on global numbering with SFC ordering */
  for (i = 0; i < n_cells; i++)
    cell_rank[i] = ((cell_num[i] - 1) / block_size) * rank_step;

  cell_io_num = fvm_io_num_destroy(cell_io_num);
}

/*----------------------------------------------------------------------------
 * Compute default face destination rank array in case of isolated faces.
 *
 * parameters:
 *   mr           <-> pointer to mesh reader helper structure
 *   comm         <-- associated MPI communicator
 *
 * returns:
 *  default rank array for faces (>= for isolated faces)
 *----------------------------------------------------------------------------*/

static int *
_default_face_rank(_mesh_reader_t     *mr,
                   MPI_Comm            comm)
{
  cs_lnum_t i;
  cs_block_dist_info_t free_face_bi;

  int n_ranks = 0, rank_id = -1;

  cs_lnum_t _n_faces = 0, n_free_faces = 0;
  cs_gnum_t _n_g_free_faces = 0, n_g_free_faces = 0;

  cs_lnum_t *free_face_ids = NULL;
  cs_coord_t *free_face_centers = NULL;

  fvm_io_num_t *free_face_io_num = NULL;
  const cs_gnum_t *free_face_num = NULL;

  int *default_rank = NULL;

  /* Initialization */

  assert((sizeof(cs_lnum_t) == 4) || (sizeof(cs_lnum_t) == 8));

  /* Count number of isolated faces */

  _n_faces = mr->face_bi.gnum_range[1] - mr->face_bi.gnum_range[0];
  n_free_faces = 0;

  for (i = 0; i < _n_faces; i++) {
    if (mr->face_cells[i*2] == 0 && mr->face_cells[i*2 + 1] == 0)
      n_free_faces += 1;
  }

  _n_g_free_faces = n_free_faces;
  MPI_Allreduce(&_n_g_free_faces, &n_g_free_faces, 1,
                CS_MPI_GNUM, MPI_SUM, comm);

  /* Return if we do not have isolated faces */

  if (n_g_free_faces == 0)
    return NULL;

  /* Initialize rank info */

  MPI_Comm_size(comm, &n_ranks);
  MPI_Comm_size(comm, &rank_id);
  free_face_bi = cs_block_dist_compute_sizes(rank_id,
                                             n_ranks,
                                             0,
                                             0,
                                             n_g_free_faces);

  /* Define distribution of isolated faces based on sfc */

  BFT_MALLOC(default_rank, _n_faces, int);
  for (i = 0; i < _n_faces; i++)
    default_rank[i] = -1;

  BFT_MALLOC(free_face_ids, n_free_faces, cs_lnum_t);
  BFT_MALLOC(free_face_centers, n_free_faces*3, cs_coord_t);

  n_free_faces = 0;
  for (i = 0; i < _n_faces; i++) {
    if (mr->face_cells[i*2] == 0 && mr->face_cells[i*2 + 1] == 0)
      free_face_ids[n_free_faces++] = i;
  }

  _precompute_free_face_center(mr,
                               n_free_faces,
                               free_face_ids,
                               free_face_centers,
                               comm);

  free_face_io_num = fvm_io_num_create_from_sfc(free_face_centers,
                                                3,
                                                n_free_faces,
                                                _sfc_type);

  BFT_FREE(free_face_centers);

  free_face_num = fvm_io_num_get_global_num(free_face_io_num);

  /* Determine rank based on global numbering with SFC ordering */
  for (i = 0; i < n_free_faces; i++) {
    default_rank[free_face_ids[i]]
      =    ((free_face_num[i] - 1) / free_face_bi.block_size)
         * free_face_bi.rank_step;
  }

  free_face_io_num = fvm_io_num_destroy(free_face_io_num);
  BFT_FREE(free_face_ids);

  return default_rank;
}

/*----------------------------------------------------------------------------
 * Organize data read by blocks in parallel and build most mesh structures.
 *
 * parameters:
 *   mesh         <-> pointer to mesh structure
 *   mesh_builder <-> pointer to mesh builder structure
 *   mr           <-> pointer to mesh reader helper structure
 *   comm         <-- associated MPI communicator
 *----------------------------------------------------------------------------*/

static void
_decompose_data_g(cs_mesh_t          *mesh,
                  cs_mesh_builder_t  *mesh_builder,
                  _mesh_reader_t     *mr,
                  MPI_Comm            comm)
{
  cs_lnum_t i;
  int n_ranks = 0;

  cs_datatype_t lnum_type = (sizeof(cs_lnum_t) == 8) ? CS_INT64 : CS_INT32;
  cs_datatype_t gnum_type = (sizeof(cs_gnum_t) == 8) ? CS_UINT64 : CS_UINT32;
  cs_datatype_t real_type = (sizeof(cs_real_t) == 8) ? CS_DOUBLE : CS_FLOAT;

  int use_cell_rank = 0;

  cs_lnum_t _n_faces = 0;
  cs_gnum_t *_face_num = NULL;
  cs_gnum_t *_face_gcells = NULL;
  cs_gnum_t *_face_gvertices = NULL;

  cs_lnum_t *_face_cells = NULL;
  cs_lnum_t *_face_gc_id = NULL;
  cs_lnum_t *_face_vertices_idx = NULL;
  cs_lnum_t *_face_vertices = NULL;

  int  *default_face_rank = NULL;
  char *face_type = NULL;
  cs_interface_set_t *face_ifs = NULL;

  cs_block_to_part_t *d = NULL;

  /* Initialization */

  MPI_Comm_size(comm, &n_ranks);

  assert((sizeof(cs_lnum_t) == 4) || (sizeof(cs_lnum_t) == 8));

  /* Different handling of cells depending on whether decomposition
     data is available or not. */

  if (mr->read_cell_rank != 0)
    use_cell_rank = 1;

  else if (_use_sfc == true && mr->read_cell_rank == 0) {

    cs_lnum_t _n_cells = mr->cell_bi.gnum_range[1] - mr->cell_bi.gnum_range[0];

    BFT_MALLOC(mr->cell_rank, _n_cells, cs_lnum_t);

    _cell_rank_by_sfc(mr,  mr->cell_rank, comm);

    use_cell_rank = 1;
  }

  if (use_cell_rank != 0) {

    d = cs_block_to_part_create_by_rank(comm,
                                        mr->cell_bi,
                                        mr->cell_rank);

    mesh->n_cells = cs_block_to_part_get_n_part_ents(d);

    BFT_MALLOC(mesh->cell_family, mesh->n_cells, cs_lnum_t);

    cs_block_to_part_copy_array(d,
                                lnum_type,
                                1,
                                mr->cell_gc_id,
                                mesh->cell_family);

    BFT_FREE(mr->cell_gc_id);

    mesh->global_cell_num = cs_block_to_part_transfer_gnum(d);

    cs_block_to_part_destroy(&d);

  }
  else {

    mesh->n_cells = mr->cell_bi.gnum_range[1] - mr->cell_bi.gnum_range[0];

    BFT_MALLOC(mesh->global_cell_num, mesh->n_cells, cs_gnum_t);

    for (i = 0; i < mesh->n_cells; i++)
      mesh->global_cell_num[i] = mr->cell_bi.gnum_range[0] + i;

    mesh->cell_family = mr->cell_gc_id;
    mr->cell_gc_id = NULL;
  }

  if (mesh->n_cells == 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Number of cells on rank %d is zero.\n"
                "(number of cells / number of processes ratio too low)."),
              (int)cs_glob_rank_id);

  /* Distribute faces */
  /*------------------*/

  default_face_rank = _default_face_rank(mr, comm);

  d = cs_block_to_part_create_by_adj_s(comm,
                                       mr->face_bi,
                                       mr->cell_bi,
                                       2,
                                       mr->face_cells,
                                       mr->cell_rank,
                                       default_face_rank);

  if (default_face_rank != NULL)
    BFT_FREE(default_face_rank);

  BFT_FREE(mr->cell_rank); /* Not needed anymore */

  _n_faces = cs_block_to_part_get_n_part_ents(d);

  BFT_MALLOC(_face_gcells, _n_faces*2, cs_gnum_t);

  /* Face -> cell connectivity */

  cs_block_to_part_copy_array(d,
                              gnum_type,
                              2,
                              mr->face_cells,
                              _face_gcells);

  BFT_FREE(mr->face_cells);

  /* Now convert face -> cell connectivity to local cell numbers */

  BFT_MALLOC(_face_cells, _n_faces*2, cs_lnum_t);

  cs_block_to_part_global_to_local(_n_faces*2,
                                   1,
                                   mesh->n_cells,
                                   mesh->global_cell_num,
                                   _face_gcells,
                                   _face_cells);

  BFT_FREE(_face_gcells);

  /* Face family */

  BFT_MALLOC(_face_gc_id, _n_faces, cs_lnum_t);

  cs_block_to_part_copy_array(d,
                              lnum_type,
                              1,
                              mr->face_gc_id,
                              _face_gc_id);

  BFT_FREE(mr->face_gc_id);

  /* Face connectivity */

  BFT_MALLOC(_face_vertices_idx, _n_faces + 1, cs_lnum_t);

  cs_block_to_part_copy_index(d,
                              mr->face_vertices_idx,
                              _face_vertices_idx);

  BFT_MALLOC(_face_gvertices, _face_vertices_idx[_n_faces], cs_gnum_t);

  cs_block_to_part_copy_indexed(d,
                                gnum_type,
                                mr->face_vertices_idx,
                                mr->face_vertices,
                                _face_vertices_idx,
                                _face_gvertices);

  BFT_FREE(mr->face_vertices_idx);
  BFT_FREE(mr->face_vertices);

  _face_num = cs_block_to_part_transfer_gnum(d);

  cs_block_to_part_destroy(&d);

  /* Vertices */

  d = cs_block_to_part_create_adj(comm,
                                  mr->vertex_bi,
                                  _face_vertices_idx[_n_faces],
                                  _face_gvertices);

  mesh->n_vertices = cs_block_to_part_get_n_part_ents(d);

  BFT_MALLOC(mesh->vtx_coord, mesh->n_vertices*3, cs_real_t);

  cs_block_to_part_copy_array(d,
                              real_type,
                              3,
                              mr->vertex_coords,
                              mesh->vtx_coord);

  BFT_FREE(mr->vertex_coords);

  mesh->global_vtx_num = cs_block_to_part_transfer_gnum(d);

  cs_block_to_part_destroy(&d);

  /* Now convert face -> vertex connectivity to local vertex numbers */

  BFT_MALLOC(_face_vertices, _face_vertices_idx[_n_faces], cs_lnum_t);

  cs_block_to_part_global_to_local(_face_vertices_idx[_n_faces],
                                   1,
                                   mesh->n_vertices,
                                   mesh->global_vtx_num,
                                   _face_gvertices,
                                   _face_vertices);

  BFT_FREE(_face_gvertices);

  /* In case of periodicity, build a cs_interface so as to obtain
     periodic face correspondants in local numbering (periodic couples
     need not be defined by the ranks owning one of the 2 members
     for the interface to be built correctly). */

  face_ifs
    = cs_interface_set_create(_n_faces,
                               NULL,
                               _face_num,
                               mesh->periodicity,
                               mr->n_perio,
                               mr->periodicity_num,
                               mr->n_per_face_couples,
                               (const cs_gnum_t *const *)mr->per_face_couples);

  /* We may now separate interior from boundary faces */

  BFT_MALLOC(face_type, _n_faces, char);

  _face_type_g(mesh,
               _n_faces,
               face_ifs,
               _face_cells,
               _face_vertices_idx,
               face_type);

  _extract_face_cell(mesh, _n_faces, _face_cells, face_type);

  {
    cs_gnum_t _n_g_free_faces = mesh->n_g_free_faces;
    MPI_Allreduce(&_n_g_free_faces, &(mesh->n_g_free_faces), 1,
                  CS_MPI_GNUM, MPI_SUM, comm);
  }

  BFT_FREE(_face_cells);

  if (mr->n_perio == 0)
    cs_interface_set_destroy(&face_ifs);

  _extract_face_vertices(mesh,
                         _n_faces,
                         _face_vertices_idx,
                         _face_vertices,
                         face_type);

  BFT_FREE(_face_vertices_idx);
  BFT_FREE(_face_vertices);

  _extract_face_gnum(mesh,
                     _n_faces,
                     _face_num,
                     face_type);

  BFT_FREE(_face_num);

  if (mr->n_perio > 0) {
    _face_ifs_to_interior(face_ifs, _n_faces, face_type);
    _extract_periodic_faces_g(mesh,
                              mesh_builder,
                              face_ifs);
    cs_interface_set_destroy(&face_ifs);
  }

  _extract_face_gc_id(mesh,
                      _n_faces,
                      _face_gc_id,
                      face_type);

  BFT_FREE(_face_gc_id);

  BFT_FREE(face_type);
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Organize data read locally and build most mesh structures
 *
 * parameters:
 *   mesh         <-- pointer to mesh structure
 *   mesh_builder <-- pointer to mesh builder structure
 *   mr           <-> pointer to mesh reader helper structure
 *----------------------------------------------------------------------------*/

static void
_decompose_data_l(cs_mesh_t          *mesh,
                  cs_mesh_builder_t  *mesh_builder,
                  _mesh_reader_t     *mr)
{
  cs_lnum_t i;

  cs_lnum_t _n_faces = 0;

  cs_lnum_t *_face_cells = NULL;
  cs_lnum_t *_face_vertices_idx = NULL;
  cs_lnum_t *_face_vertices = NULL;

  char *face_type = NULL;

  /* Initialization */

  assert((sizeof(cs_lnum_t) == 4) || (sizeof(cs_lnum_t) == 8));

  mesh->n_cells = mr->cell_bi.gnum_range[1] - 1;

  /* Cell families are already of the correct type,
     so they can simply be moved */

  mesh->cell_family = mr->cell_gc_id;
  mr->cell_gc_id = NULL;

  /* Build faces */
  /*-------------*/

  _n_faces = mr->face_bi.gnum_range[1] - 1;

  /* Now copy face -> cell connectivity to local cell numbers */

  BFT_MALLOC(_face_cells, _n_faces*2, cs_lnum_t);

  for (i = 0; i < _n_faces; i++) {
    _face_cells[i*2] = mr->face_cells[i*2];
    _face_cells[i*2 + 1] = mr->face_cells[i*2 + 1];
  }

  BFT_FREE(mr->face_cells);

  /* Face connectivity */

  BFT_MALLOC(_face_vertices_idx, _n_faces + 1, cs_lnum_t);

  for (i = 0; i < _n_faces+1; i++)
    _face_vertices_idx[i] = mr->face_vertices_idx[i];

  BFT_FREE(mr->face_vertices_idx);

  BFT_MALLOC(_face_vertices, _face_vertices_idx[_n_faces], cs_lnum_t);

  for (i = 0; i < _face_vertices_idx[_n_faces]; i++)
    _face_vertices[i] = mr->face_vertices[i];

  BFT_FREE(mr->face_vertices);

  /* Vertices */

  mesh->n_vertices = mr->vertex_bi.gnum_range[1] - 1;

  mesh->vtx_coord = mr->vertex_coords;
  mr->vertex_coords = NULL;

  /* We may now separate interior from boundary faces */

  BFT_MALLOC(face_type, _n_faces, char);

  _face_type_l(mesh,
               _n_faces,
               mr->n_per_face_couples,
               (const cs_gnum_t *const *)mr->per_face_couples,
               _face_cells,
               _face_vertices_idx,
               face_type);

  _extract_face_cell(mesh, _n_faces, _face_cells, face_type);

  BFT_FREE(_face_cells);

  if (mr->n_perio > 0) {

    /* Transfer arrays from reader to builder, then renumber couples */

    _extract_periodic_faces_l(mr,
                              mesh_builder,
                              mesh->n_init_perio,
                              _n_faces,
                              face_type);
  }

  _extract_face_vertices(mesh,
                         _n_faces,
                         _face_vertices_idx,
                         _face_vertices,
                         face_type);

  BFT_FREE(_face_vertices_idx);
  BFT_FREE(_face_vertices);

  _extract_face_gc_id(mesh,
                      _n_faces,
                      mr->face_gc_id,
                      face_type);

  BFT_FREE(mr->face_gc_id);

  BFT_FREE(face_type);
}

/*----------------------------------------------------------------------------
 * Rename groups in a mesh.
 *
 * parameters:
 *   mesh            <-> mesh being modified
 *   start_id        <-- id of first group to rename
 *   n_group_renames <-- number of group rename couples
 *   old_group_names <-- old group names (size: n_group_renames)
 *   new_group_names <-- new group names (size: n_group_renames)
 *----------------------------------------------------------------------------*/

static void
_mesh_groups_rename(cs_mesh_t          *mesh,
                    size_t              start_id,
                    size_t              n_group_renames,
                    const char  *const *old_group_names,
                    const char  *const *new_group_names)
{
  size_t  i, j, k;
  int  have_rename = 0;
  size_t  end_id = mesh->n_groups;
  size_t  n_ids = 0;
  int  *new_group_id = NULL;

  if (end_id > start_id)
    n_ids = end_id - start_id;

  /* Check for rename matches */

  BFT_MALLOC(new_group_id, n_ids, int);

  for (i = 0, j = start_id; i < n_ids; i++, j++) {
    const char *g_name = mesh->group_lst + (mesh->group_idx[j] - 1);
    new_group_id[i] = -1;
    for (k = 0; k < n_group_renames; k++) {
      if (strcmp(g_name, old_group_names[k]) == 0) {
        new_group_id[i] = k;
        have_rename = 1;
        break;
      }
    }
  }

  /* Now rename matched groups */

  if (have_rename) {

    size_t new_size = 0;
    size_t old_size = mesh->group_idx[end_id] - mesh->group_idx[start_id];
    int *saved_idx = NULL;
    char *saved_names = NULL;

    BFT_MALLOC(saved_idx, n_ids + 1, int);
    BFT_MALLOC(saved_names, old_size, char);

    for (i = 0; i < n_ids+1; i++)
      saved_idx[i] = mesh->group_idx[start_id + i] - mesh->group_idx[start_id];
    memcpy(saved_names,
           mesh->group_lst + (mesh->group_idx[start_id] - 1),
           old_size);

    /* Update index */

    for (i = 0; i < n_ids; i++) {
      const char *new_src = NULL;
      new_size += 1;
      if (new_group_id[i] > -1)
        new_src = new_group_names[new_group_id[i]];
      else
        new_src = saved_names + saved_idx[i];
      if (new_src != NULL)
        new_size += strlen(new_src);
      mesh->group_idx[start_id + i + 1]
        = mesh->group_idx[start_id] + new_size;
    }

    BFT_REALLOC(mesh->group_lst, mesh->group_idx[mesh->n_groups] - 1, char);

    for (i = 0, j = start_id; i < n_ids; i++, j++) {
      char *new_dest = mesh->group_lst + (mesh->group_idx[j] - 1);
      const char *new_src = NULL;
      if (new_group_id[i] > -1)
        new_src = new_group_names[new_group_id[i]];
      else
        new_src = saved_names + saved_idx[i];
      if (new_src != NULL)
        strcpy(new_dest, new_src);
      else
        strcpy(new_dest, "");
    }

    BFT_FREE(saved_names);
    BFT_FREE(saved_idx);

    /* Set mesh modification flag */

    mesh->modified = 1;

  }

  BFT_FREE(new_group_id);
}

/*----------------------------------------------------------------------------
 * Read sections from the pre-processor about the dimensions of mesh
 *
 * This function updates the information in the mesh and mesh reader
 * structures relative to the data in the given file.
 *
 * parameters
 *   mesh <-> pointer to mesh structure
 *   n_gc <-- number of group classes last read
 *----------------------------------------------------------------------------*/

static void
_colors_to_groups(cs_mesh_t  *mesh,
                  int         n_gc)
{
  cs_int_t  i, j;
  int  n_colors = 0;
  int  color_names_size = 0;

  /* Counting pass */

  for (j = 0; j < mesh->n_max_family_items; j++) {
    for (i = mesh->n_families - n_gc; i < mesh->n_families; i++) {
      if (mesh->family_item[mesh->n_families*j + i] > 0) {
        int color_id = mesh->family_item[mesh->n_families*j + i];
        int name_size = 1;
        while (color_id > 0) {
          color_id /= 10;
          name_size += 1;
        }
        n_colors += 1;
        color_names_size += name_size;
      }
    }
  }

  /* Reallocation */

  if (n_colors > 0) {
    if (mesh->n_groups > 0) {
      BFT_REALLOC(mesh->group_idx, mesh->n_groups + n_colors + 1, cs_int_t);
      BFT_REALLOC(mesh->group_lst,
                  mesh->group_idx[mesh->n_groups] - 1 + color_names_size,
                  char);
    }
    else {
      BFT_MALLOC(mesh->group_idx, n_colors + 1, cs_int_t);
      BFT_MALLOC(mesh->group_lst, color_names_size, char);
      mesh->group_idx[0] = 1;
    }
  }

  /* Assignment */

  for (j = 0; j < mesh->n_max_family_items; j++) {
    for (i = mesh->n_families - n_gc; i < mesh->n_families; i++) {
      if (mesh->family_item[mesh->n_families*j + i] > 0) {
        int color_id = mesh->family_item[mesh->n_families*j + i];
        int name_size = 1;
        cs_int_t group_lst_end = mesh->group_idx[mesh->n_groups] - 1;
        sprintf(mesh->group_lst + group_lst_end, "%d", color_id);
        while (color_id > 0) {
          color_id /= 10;
          name_size += 1;
        }
        mesh->group_idx[mesh->n_groups + 1]
          = mesh->group_idx[mesh->n_groups] + name_size;
        mesh->n_groups += 1;
        mesh->family_item[mesh->n_families*j + i] = - mesh->n_groups;
      }
    }
  }

}

/*----------------------------------------------------------------------------
 * Read sections from the pre-processor about the dimensions of mesh
 *
 * This function updates the information in the mesh and mesh reader
 * structures relative to the data in the given file.
 *
 * parameters
 *   mesh     <-> pointer to mesh structure
 *   mr       <-> pointer to mesh reader structure
 *   file_id  <-- file id in mesh reader
 *----------------------------------------------------------------------------*/

static void
_read_dimensions(cs_mesh_t       *mesh,
                 _mesh_reader_t  *mr,
                 int              file_id)
{
  cs_int_t  i, j;
  cs_io_sec_header_t  header;

  cs_gnum_t n_elts = 0;
  int        n_gc = 0;
  int        n_gc_props_max = 0;
  int        n_groups = 0;
  int        n_init_perio = 0;
  bool       dim_read = false;
  bool       end_read = false;
  cs_io_t   *pp_in = NULL;

  _mesh_file_info_t  *f = NULL;

  const char  *unexpected_msg = N_("Section of type <%s> on <%s>\n"
                                   "unexpected or of incorrect size");

  if (file_id < 0 || file_id >= mr->n_files)
    return;

  f = mr->file_info + file_id;
  f->offset = 0;

  mr->gc_id_shift[file_id] = mesh->n_families;

  /* Initialize reading of Preprocessor output */

  bft_printf(_(" Reading metadata from file: \"%s\"\n"), f->filename);

#if defined(HAVE_MPI)
  pp_in = cs_io_initialize(f->filename,
                           "Face-based mesh definition, R0",
                           CS_IO_MODE_READ,
                           CS_FILE_NO_MPI_IO,
                           CS_IO_ECHO_NONE,
                           cs_glob_mpi_comm);
#else
  pp_in = cs_io_initialize(f->filename,
                           "Face-based mesh definition, R0",
                           CS_IO_MODE_READ,
                           CS_IO_ECHO_NONE,
                           -1);
#endif

  /* Loop on read sections */

  while (end_read == false) {

    /* Receive headers and clean header names */

    cs_io_read_header(pp_in, &header);

    /* Treatment according to the header name */

    if (strncmp(header.sec_name, "EOF", CS_IO_NAME_LEN)
        == 0) {
      cs_io_finalize(&pp_in);
      pp_in = NULL;
    }

    if (strncmp(header.sec_name, "start_block:dimensions",
                CS_IO_NAME_LEN) == 0) {

      if (dim_read == false)
        dim_read = true;
      else
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));

    }
    else if (strncmp(header.sec_name, "end_block:dimensions",
                     CS_IO_NAME_LEN) == 0) {

      if (dim_read == true) {
        dim_read = false;
        end_read = true;
      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));

    }

    /* Receive dimensions from the pre-processor */

    else if (strncmp(header.sec_name, "n_cells",
                     CS_IO_NAME_LEN) == 0) {

      if (dim_read != true || header.n_vals != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else {
        cs_gnum_t _n_g_cells;
        cs_io_set_cs_gnum(&header, pp_in);
        cs_io_read_global(&header, &_n_g_cells, pp_in);
        mesh->n_g_cells += _n_g_cells;
      }

    }
    else if (strncmp(header.sec_name, "n_faces",
                     CS_IO_NAME_LEN) == 0) {

      if (dim_read != true || header.n_vals != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else {
        cs_gnum_t _n_g_faces;
        cs_io_set_cs_gnum(&header, pp_in);
        cs_io_read_global(&header, &_n_g_faces, pp_in);
        mr->n_g_faces += _n_g_faces;
      }

    }
    else if (strncmp(header.sec_name, "n_vertices",
                     CS_IO_NAME_LEN) == 0) {

      if (dim_read != true || header.n_vals != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else {
        cs_gnum_t _n_g_vertices;
        cs_io_set_cs_gnum(&header, pp_in);
        cs_io_read_global(&header, &_n_g_vertices, pp_in);
        mesh->n_g_vertices += _n_g_vertices;
      }

    }
    else if (strncmp(header.sec_name, "face_vertices_size",
                     CS_IO_NAME_LEN) == 0) {

      if (dim_read != true || header.n_vals != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else {
        cs_gnum_t _n_g_face_connect_size;
        cs_io_set_cs_gnum(&header, pp_in);
        cs_io_read_global(&header, &_n_g_face_connect_size, pp_in);
        mr->n_g_face_connect_size += _n_g_face_connect_size;
      }

    }
    else if (strncmp(header.sec_name, "n_group_classes",
                     CS_IO_NAME_LEN) == 0) {

      if (dim_read != true || header.n_vals != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else {
        cs_io_set_cs_lnum(&header, pp_in);
        cs_io_read_global(&header, &n_gc, pp_in);
        mesh->n_families += n_gc;
      }

    }
    else if (strncmp(header.sec_name, "n_group_class_props_max",
                     CS_IO_NAME_LEN) == 0) {

      if (dim_read != true || header.n_vals != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else {
        cs_io_set_cs_lnum(&header, pp_in);
        cs_io_read_global(&header, &n_gc_props_max, pp_in);
        if (n_gc_props_max > mesh->n_max_family_items) {
          /* Update (pad) previous definitions */
          BFT_REALLOC(mesh->family_item,
                      mesh->n_families*n_gc_props_max,
                      cs_int_t);
          for (i = mesh->n_max_family_items;
               i < n_gc_props_max;
               i++) {
            for (j = 0; j < mesh->n_families - n_gc; j++)
              mesh->family_item[(mesh->n_families - n_gc)*i + j] = 0;
          }
          mesh->n_max_family_items = n_gc_props_max;
        }
      }

    }
    else if (strncmp(header.sec_name, "n_groups",
                     CS_IO_NAME_LEN) == 0) {

      if (dim_read != true || header.n_vals != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else {
        cs_io_set_cs_lnum(&header, pp_in);
        cs_io_read_global(&header, &n_groups, pp_in);
        mesh->n_groups += n_groups;
      }

    }
    else if (strncmp(header.sec_name, "group_name_index",
                     CS_IO_NAME_LEN) == 0) {

      if ((cs_int_t)header.n_vals != n_groups + 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else {
        cs_io_set_cs_lnum(&header, pp_in);
        if (mesh->group_idx == NULL) {
          BFT_MALLOC(mesh->group_idx, mesh->n_groups + 1, cs_int_t);
          cs_io_read_global(&header, mesh->group_idx, pp_in);
        }
        else {
          cs_int_t *_group_idx = NULL;
          BFT_REALLOC(mesh->group_idx, mesh->n_groups + 1, cs_int_t);
          BFT_MALLOC(_group_idx, n_groups + 1, cs_int_t);
          cs_io_read_global(&header, _group_idx, pp_in);
          for (i = 0, j = mesh->n_groups - n_groups; i < n_groups; i++, j++)
            mesh->group_idx[j + 1] = (   mesh->group_idx[j]
                                      + _group_idx[i+1] - _group_idx[i]);
          BFT_FREE(_group_idx);
        }
      }

    }
    else if (strncmp(header.sec_name, "group_name",
                     CS_IO_NAME_LEN) == 0) {

      if (   mesh->group_idx == NULL
          || (  (cs_int_t)header.n_vals
              != (  mesh->group_idx[mesh->n_groups]
                  - mesh->group_idx[mesh->n_groups - n_groups])))
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else {
        i = mesh->group_idx[mesh->n_groups - n_groups] - mesh->group_idx[0];
        BFT_REALLOC(mesh->group_lst, i + header.n_vals + 1, char);
        cs_io_read_global(&header, mesh->group_lst + i, pp_in);
      }

    }
    else if (strncmp(header.sec_name, "group_class_properties",
                     CS_IO_NAME_LEN) == 0) {

      n_elts = mesh->n_families * mesh->n_max_family_items;
      if (dim_read != true || header.n_vals != n_gc*n_gc_props_max)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else {
        if (mesh->family_item == NULL)
          BFT_MALLOC(mesh->family_item, n_elts, cs_int_t);
        if (mesh->n_families == n_gc)
          cs_io_read_global(&header, mesh->family_item, pp_in);
        else {
          cs_int_t *_family_item = NULL;
          BFT_REALLOC(mesh->family_item, n_elts, cs_int_t);
          BFT_MALLOC(_family_item, header.n_vals, cs_int_t);
          cs_io_read_global(&header, _family_item, pp_in);
          /* Shift previous data */
          for (j = mesh->n_max_family_items - 1; j > 0; j--) {
            for (i = mesh->n_families - n_gc - 1; i > -1; i--)
              mesh->family_item[mesh->n_families*j + i]
                = mesh->family_item[(mesh->n_families - n_gc)*j + i];
          }
          for (i = 0; i < n_gc; i++) {
            /* Copy group class data, shifting group names if necessary */
            for (j = 0; j < n_gc_props_max; j++) {
              int _family_item_j = _family_item[n_gc*j + i];
              if (_family_item_j < 0)
                _family_item_j -= (mesh->n_groups - n_groups);
              mesh->family_item[mesh->n_families*j + (mesh->n_families-n_gc+i)]
                = _family_item_j;
            }
            /* Pad if necessary */
            for (j = n_gc_props_max; j < mesh->n_max_family_items; j++)
              mesh->family_item[mesh->n_families*j + (mesh->n_families-n_gc+i)]
                = 0;
          }
          BFT_FREE(_family_item);
        }
        /* Transform colors to group names if present */
        _colors_to_groups(mesh, n_gc);
      }

      if (f->n_group_renames > 0)
        _mesh_groups_rename(mesh,
                            mesh->n_groups - n_groups,
                            f->n_group_renames,
                            f->old_group_names,
                            f->new_group_names);

    }

    /* Additional sections for periodicity. */

    else if (strncmp(header.sec_name, "n_periodic_directions",
                     CS_IO_NAME_LEN) == 0) {

      if (dim_read != true || header.n_vals != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else {

        cs_io_set_cs_lnum(&header, pp_in);
        cs_io_read_global(&header, &n_init_perio, pp_in);
        mesh->n_init_perio += n_init_perio;

        assert(n_init_perio > 0);

        if (mesh->periodicity == NULL)
          mesh->periodicity = fvm_periodicity_create(0.001);

        BFT_REALLOC(mr->periodicity_num, mesh->n_init_perio, int);
        BFT_REALLOC(mr->n_per_face_couples, mesh->n_init_perio, cs_lnum_t);
        BFT_REALLOC(mr->n_g_per_face_couples, mesh->n_init_perio, cs_gnum_t);
        BFT_REALLOC(mr->per_face_couples, mesh->n_init_perio, cs_gnum_t *);

        mr->n_perio = mesh->n_init_perio;

        for (i = mesh->n_init_perio - n_init_perio;
             i < mesh->n_init_perio;
             i++) {
          mr->periodicity_num[i] = i+1;
          mr->n_per_face_couples[i] = 0;
          mr->n_g_per_face_couples[i] = 0;
          mr->per_face_couples[i] = NULL;
        }
      }

    }
    else if (strncmp(header.sec_name, "n_periodic_rotations",
                     CS_IO_NAME_LEN) == 0) {

      if (dim_read != true || header.n_vals != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else {
        cs_int_t n_rot_perio;
        cs_io_set_cs_lnum(&header, pp_in);
        cs_io_read_global(&header, &n_rot_perio, pp_in);
        if (n_rot_perio > 0)
          mesh->have_rotation_perio = 1;
      }

    }
    else if (strncmp(header.sec_name, "n_periodic_faces",
                     CS_IO_NAME_LEN) == 0) {

      if ((cs_int_t)header.n_vals != n_init_perio)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else {
        size_t dest_offset = mesh->n_init_perio - n_init_perio;
        cs_io_set_cs_gnum(&header, pp_in);
        cs_io_read_global(&header,
                          mr->n_g_per_face_couples + dest_offset,
                          pp_in);
        for (i = dest_offset; i < mesh->n_init_perio; i++)
          mr->n_g_per_face_couples[i] /= 2;
      }

    }
    else
      bft_error(__FILE__, __LINE__, 0,
                _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));

  } /* End of test on headers */

  /* Close file */

  f->offset = cs_io_get_offset(pp_in);
  cs_io_finalize(&pp_in);
}

/*----------------------------------------------------------------------------
 * Compute value range information for a given data section.
 *
 * parameters:
 *   header          <-- pointer to current file section header data
 *   pp_in           <-- pointer to current file
 *   n_g_elts        <-> global number of elements
 *   n_g_elts_read   <-> global number of elements already read
 *   n_location_vals <-- number of values for each location
 *   is_index        <-- 1 if data is an index, 0 otherwise
 *   gnum_range      <-- global number range for all elements on this rank
 *   gnum_range_cur  --> global number range for elements in current file
 *   n_g_elts_cur    --> global number of elements in current file
 *   n_vals          --> expected number of local values from all files
 *   n_vals_cur      --> number of local values from current file
 *----------------------------------------------------------------------------*/

static void
_data_range(cs_io_sec_header_t  *header,
            const cs_io_t       *pp_in,
            cs_gnum_t            n_g_elts,
            cs_gnum_t            n_g_elts_read,
            size_t               n_location_vals,
            size_t               is_index,
            const cs_gnum_t      gnum_range[2],
            cs_gnum_t            gnum_range_cur[2],
            cs_gnum_t           *n_g_elts_cur,
            cs_lnum_t           *n_vals,
            cs_lnum_t           *n_vals_cur)
{
  size_t i;

  /* Initialization */

  gnum_range_cur[0] = gnum_range[0];
  gnum_range_cur[1] = gnum_range[1];

  *n_g_elts_cur = (header->n_vals - is_index) / n_location_vals;
  *n_vals = (gnum_range[1] - gnum_range[0]) * n_location_vals;
  *n_vals_cur = 0;

  if (*n_g_elts_cur + n_g_elts_read > n_g_elts)
    bft_error(__FILE__, __LINE__, 0,
              _("Section of type <%s> on <%s>\n"
                "has incorrect size (current: %llu, read: %llu, total: %llu."),
              header->sec_name, cs_io_get_name(pp_in),
              (unsigned long long)(*n_g_elts_cur),
              (unsigned long long)n_g_elts_read,
              (unsigned long long)n_g_elts);

  else if (header->n_location_vals != n_location_vals)
    bft_error(__FILE__, __LINE__, 0,
              _("Section of type <%s> on <%s>\n"
                "has incorrect number of values per location."),
              header->sec_name, cs_io_get_name(pp_in));

  else {

    /* Compute range for this file (based on range for all files,
       and parts of this range already read) */

    for (i = 0; i < 2; i++) {
      if (gnum_range_cur[i] <= n_g_elts_read)
        gnum_range_cur[i] = 1;
      else
        gnum_range_cur[i] -= n_g_elts_read;
      if (gnum_range_cur[i] > *n_g_elts_cur)
        gnum_range_cur[i] = *n_g_elts_cur + 1;
    }

    if (gnum_range[1] > gnum_range[0])
      *n_vals_cur = (gnum_range_cur[1] - gnum_range_cur[0]) * n_location_vals;
  }

  /* Index adds past-the-last value */
  if (is_index == 1) {
    *n_vals += 1;
    *n_vals_cur += 1;
  }
}

/*----------------------------------------------------------------------------
 * Transform coordinates.
 *
 * parameters:
 *   n_coords <-- number of coordinates.
 *   coords   <-> array of coordinates
 *   matrix   <-- matrix of the transformation in homogeneous coord.
 *                last line = [0; 0; 0; 1] (Not used here)
 *----------------------------------------------------------------------------*/

static void
_transform_coords(size_t         n_coords,
                  cs_real_t     *coords,
                  const double   matrix[])

{
  size_t  i;
  double  _matrix[3][4];
  memcpy(_matrix, matrix, 12*sizeof(double));

  for (i = 0; i < n_coords; i++) {

    size_t j, coo_shift;
    cs_real_t  tmp_coord[3];

    coo_shift = i*3;

    for (j = 0; j < 3; j++)
      tmp_coord[j] = coords[coo_shift + j];

    for (j = 0; j < 3; j++) {
      coords[coo_shift + j] = (  _matrix[j][0] * tmp_coord[0]
                               + _matrix[j][1] * tmp_coord[1]
                               + _matrix[j][2] * tmp_coord[2]
                               + _matrix[j][3]);
    }

  }
}

/*----------------------------------------------------------------------------
 * Invert a homogeneous transformation matrix.
 *
 * parameters:
 *   a  <-- matrix
 *   b  --> inverse matrix
 *----------------------------------------------------------------------------*/

static void
_inverse_transf_matrix(double  a[3][4],
                       double  b[3][4])
{
  int i, j, k, k_pivot;

  double abs_pivot, abs_a_ki, factor;
  double _a[3][4];
  double _a_swap[3], _b_swap[3];

  /* Copy array (avoid modifying a), and initialize inverse */

  memcpy(_a, a, sizeof(double)*12);

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 4; j++)
      b[i][j] = 0.0;
  }
  for (i = 0; i < 3; i++)
    b[i][i] = 1.0;

  /* Forward elimination */

  for (i = 0; i < 2; i++) {

    /* Seek best pivot */

    k_pivot = i;
    abs_pivot = fabs(_a[i][i]);

    for (k = i+1; k < 3; k++) {
      abs_a_ki = fabs(_a[k][i]);
      if (abs_a_ki > abs_pivot) {
        abs_pivot = abs_a_ki;
        k_pivot = k;
      }
    }

    /* Swap if necessary */

    if (k_pivot != i) {
      for (j = 0; j < 4; j++) {
        _a_swap[j] = _a[i][j];
        _a[i][j] = _a[k_pivot][j];
        _a[k_pivot][j] = _a_swap[j];
        _b_swap[j] = b[i][j];
        b[i][j] = b[k_pivot][j];
        b[k_pivot][j] = _b_swap[j];
      }
    }

    if (abs_pivot < 1.0e-18)
      bft_error(__FILE__, __LINE__, 0,
                _("User-defined Transformation matrix is not invertible:\n"
                  "  [[%12.5f %12.5f %12.5f %12.5f]\n"
                  "   [%12.5f %12.5f %12.5f %12.5f]\n"
                  "   [%12.5f %12.5f %12.5f %12.5f]\n"
                  "   [%12.5f %12.5f %12.5f %12.5f]]\n"),
                a[0][0], a[0][1], a[0][2], a[0][3],
                a[1][0], a[1][1], a[1][2], a[1][3],
                a[2][0], a[2][1], a[2][2], a[2][3],
                0., 0., 0., 1.);

    /* Eliminate values */

    for (k = i+1; k < 3; k++) {
      factor = _a[k][i] / _a[i][i];
      _a[k][i] = 0.0;
      for (j = i+1; j < 4; j++) {
        _a[k][j] -= _a[i][j]*factor;
        b[k][j] -= b[i][j]*factor;
      }
    }

  }

  /* Normalize lines */

  for (i = 0; i < 3; i++) {
    factor = _a[i][i];
    for (j = 0; j < 4; j++) {
      _a[i][j] /= factor;
      b[i][j] /= factor;
    }
  }

  /* Eliminate last column */

  for (j = 3; j > 0; j--) {
    for (i = 0; i < j; i++) {
      b[i][j] -= _a[i][j];
      _a[i][j] = 0.0;
    }
  }
}

/*----------------------------------------------------------------------------
 * Combine transformation matrixes (c = a.b)
 *
 * parameters:
 *   a <-- first transformation matrix
 *   b <-- second transformation matrix
 *   c --> combined transformation matrix
 *---------------------------------------------------------------------------*/

static void
_combine_tr_matrixes(double  a[3][4],
                     double  b[3][4],
                     double  c[3][4])
{
  c[0][0] = a[0][0]*b[0][0] + a[0][1]*b[1][0] + a[0][2]*b[2][0];
  c[0][1] = a[0][0]*b[0][1] + a[0][1]*b[1][1] + a[0][2]*b[2][1];
  c[0][2] = a[0][0]*b[0][2] + a[0][1]*b[1][2] + a[0][2]*b[2][2];
  c[0][3] = a[0][0]*b[0][3] + a[0][1]*b[1][3] + a[0][2]*b[2][3] + a[0][3];

  c[1][0] = a[1][0]*b[0][0] + a[1][1]*b[1][0] + a[1][2]*b[2][0];
  c[1][1] = a[1][0]*b[0][1] + a[1][1]*b[1][1] + a[1][2]*b[2][1];
  c[1][2] = a[1][0]*b[0][2] + a[1][1]*b[1][2] + a[1][2]*b[2][2];
  c[1][3] = a[1][0]*b[0][3] + a[1][1]*b[1][3] + a[1][2]*b[2][3] + a[1][3];

  c[2][0] = a[2][0]*b[0][0] + a[2][1]*b[1][0] + a[2][2]*b[2][0];
  c[2][1] = a[2][0]*b[0][1] + a[2][1]*b[1][1] + a[2][2]*b[2][1];
  c[2][2] = a[2][0]*b[0][2] + a[2][1]*b[1][2] + a[2][2]*b[2][2];
  c[2][3] = a[2][0]*b[0][3] + a[2][1]*b[1][3] + a[2][2]*b[2][3] + a[2][3];
}

/*----------------------------------------------------------------------------
 * Transform a periodicity matrix in case of coordinates transformation
 *
 * The new transformation is the combination of the following:
 *  - switch to initial coordinates
 *  - apply initial periodicity
 *  - switch to new coordinates
 *
 * parameters:
 *   matrix       <-- matrix of the transformation in homogeneous coord.
 *                    last line = [0; 0; 0; 1] (Not used here)
 *   perio_matrix <-- matrix of the periodic matrix transformation
 *                    in homogeneous coord.
 *----------------------------------------------------------------------------*/

static void
_transform_perio_matrix(const double  matrix[],
                        double        perio_matrix[3][4])

{
  double  _m[3][4];
  double  _inv_m[3][4], _tmp_m[3][4];

  memcpy(_m, matrix, 12*sizeof(double));

  _inverse_transf_matrix(_m, _inv_m);

  _combine_tr_matrixes(perio_matrix, _inv_m, _tmp_m);
  _combine_tr_matrixes(_m, _tmp_m, perio_matrix);
}

/*----------------------------------------------------------------------------
 * Read pre-processor mesh data for a given mesh and finalize input.
 *
 * parameters:
 *   file_id <-- id of file handled by mesh builder
 *   mesh    <-> pointer to mesh structure
 *   mr      <-> pointer to mesh reader structure
 *   echo    <-- echo (verbosity) level
 *----------------------------------------------------------------------------*/

static void
_read_data(int              file_id,
           cs_mesh_t       *mesh,
           _mesh_reader_t  *mr,
           long             echo)
{
  cs_int_t  perio_id, perio_type;
  cs_io_sec_header_t  header;

  cs_real_t  perio_matrix[3][4];

  cs_int_t  perio_num = -1;
  bool       end_read = false;
  bool       data_read = false;
  cs_io_t  *pp_in = NULL;

  int gc_id_shift = mr->gc_id_shift[file_id];

  int n_perio_read = 0;
  cs_lnum_t n_cells = 0;
  cs_lnum_t n_faces = 0;
  cs_lnum_t n_vertices = 0;
  cs_lnum_t n_face_connect_size = 0;
  cs_gnum_t n_g_cells = 0;
  cs_gnum_t n_g_faces = 0;
  cs_gnum_t n_g_vertices = 0;
  cs_gnum_t n_g_face_connect_size = 0;

  cs_gnum_t face_vtx_range[2] = {0, 0};
  _mesh_file_info_t  *f = NULL;

  const char  *unexpected_msg = N_("Section of type <%s> on <%s>\n"
                                   "unexpected or of incorrect size.");

  f = mr->file_info + file_id;

#if defined(HAVE_MPI)
  pp_in = cs_io_initialize(f->filename,
                           "Face-based mesh definition, R0",
                           CS_IO_MODE_READ,
                           cs_glob_io_hints,
                           echo,
                           cs_glob_mpi_comm);
#else
  pp_in = cs_io_initialize(f->filename,
                           "Face-based mesh definition, R0",
                           CS_IO_MODE_READ,
                           echo,
                           -1);
#endif

  cs_io_set_offset(pp_in, f->offset);

  echo = cs_io_get_echo(pp_in);

  /* Loop on sections read */

  while (end_read == false) {

    /* Receive header and clean header name */

    cs_io_read_header(pp_in, &header);

    /* Process according to the header name */

    if (strncmp(header.sec_name, "EOF", CS_IO_NAME_LEN)
        == 0) {
      cs_io_finalize(&pp_in);
      pp_in = NULL;
    }

    if (strncmp(header.sec_name, "start_block:data",
                CS_IO_NAME_LEN) == 0) {

      if (data_read == false)
        data_read = true;
      else
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));

    }
    else if (strncmp(header.sec_name, "end_block:data",
                     CS_IO_NAME_LEN) == 0) {

      if (data_read == true) {
        data_read = false;
        end_read = true;
      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));

    }

    /* Read data from the pre-processor output file */

    else {

      cs_gnum_t gnum_range_cur[2];

      cs_lnum_t n_vals = 0;
      cs_lnum_t n_vals_cur = 0;
      cs_lnum_t val_offset_cur = 0;

      if (data_read != true)
        bft_error(__FILE__, __LINE__, 0,
                  _("Section of type <%s> on <%s>\n"
                    "unexpected."), header.sec_name, cs_io_get_name(pp_in));

      /* Face-cells connectivity */

      if (strncmp(header.sec_name, "face_cells", CS_IO_NAME_LEN) == 0) {

        /* Compute range for current file  */
        _data_range(&header,
                    pp_in,
                    mr->n_g_faces,
                    mr->n_g_faces_read,
                    2,
                    0,
                    mr->face_bi.gnum_range,
                    gnum_range_cur,
                    &n_g_faces,
                    &n_vals,
                    &n_vals_cur);

        n_faces = n_vals_cur / 2;
        val_offset_cur = mr->n_faces_read * 2;

        /* Allocate for first file read */
        if (mr->face_cells == NULL)
          BFT_MALLOC(mr->face_cells, n_vals, cs_gnum_t);

        /* Read data */
        cs_io_set_cs_gnum(&header, pp_in);
        cs_io_read_block(&header, gnum_range_cur[0], gnum_range_cur[1],
                         mr->face_cells + val_offset_cur, pp_in);

        /* Shift referenced cell numbers in case of appended data */
        if (mr->n_g_cells_read > 0) {
          cs_lnum_t ii;
          for (ii = 0; ii < n_vals_cur; ii++) {
            if (mr->face_cells[val_offset_cur + ii] != 0)
              mr->face_cells[val_offset_cur + ii] += mr->n_g_cells_read;
          }
        }
      }

      /* Cell group class values */

      else if (strncmp(header.sec_name, "cell_group_class_id",
                     CS_IO_NAME_LEN) == 0) {

        /* Compute range for current file  */
        _data_range(&header,
                    pp_in,
                    mesh->n_g_cells,
                    mr->n_g_cells_read,
                    1,
                    0,
                    mr->cell_bi.gnum_range,
                    gnum_range_cur,
                    &n_g_cells,
                    &n_vals,
                    &n_vals_cur);

        n_cells = n_vals_cur;
        val_offset_cur = mr->n_cells_read;

        /* Allocate for first file read */
        if (mr->cell_gc_id == NULL)
          BFT_MALLOC(mr->cell_gc_id, n_vals, cs_int_t);

        /* Read data */
        cs_io_set_cs_lnum(&header, pp_in);
        cs_io_read_block(&header, gnum_range_cur[0], gnum_range_cur[1],
                         mr->cell_gc_id + val_offset_cur, pp_in);

        /* Shift referenced numbers in case of appended data */
        if (gc_id_shift > 0) {
          cs_lnum_t ii;
          for (ii = 0; ii < n_vals_cur; ii++) {
            if (mr->cell_gc_id[val_offset_cur + ii] != 0)
              mr->cell_gc_id[val_offset_cur + ii] += gc_id_shift;
          }
        }
      }

      /* Face group class values */

      else if (strncmp(header.sec_name, "face_group_class_id",
                       CS_IO_NAME_LEN) == 0) {

        /* Compute range for current file  */
        _data_range(&header,
                    pp_in,
                    mr->n_g_faces,
                    mr->n_g_faces_read,
                    1,
                    0,
                    mr->face_bi.gnum_range,
                    gnum_range_cur,
                    &n_g_faces,
                    &n_vals,
                    &n_vals_cur);

        n_faces = n_vals_cur;
        val_offset_cur = mr->n_faces_read;

        /* Allocate for first file read */
        if (mr->face_gc_id == NULL)
          BFT_MALLOC(mr->face_gc_id, n_vals, cs_int_t);

        /* Read data */
        cs_io_set_cs_lnum(&header, pp_in);
        cs_io_read_block(&header, gnum_range_cur[0], gnum_range_cur[1],
                         mr->face_gc_id + val_offset_cur, pp_in);

        /* Shift referenced numbers in case of appended data */
        if (gc_id_shift > 0) {
          cs_lnum_t ii;
          for (ii = 0; ii < n_vals_cur; ii++) {
            if (mr->face_gc_id[val_offset_cur + ii] != 0)
              mr->face_gc_id[val_offset_cur + ii] += gc_id_shift;
          }
        }
      }

      /* Face -> vertices connectivity */

      else if (strncmp(header.sec_name, "face_vertices_index",
                       CS_IO_NAME_LEN) == 0) {

        cs_lnum_t ii;
        cs_lnum_t idx_offset_shift = 0;
        cs_gnum_t idx_gnum_shift = 0;
        cs_gnum_t *_g_face_vertices_idx = NULL;

        /* Compute range for current file  */
        _data_range(&header,
                    pp_in,
                    mr->n_g_faces,
                    mr->n_g_faces_read,
                    1,
                    1,
                    mr->face_bi.gnum_range,
                    gnum_range_cur,
                    &n_g_faces,
                    &n_vals,
                    &n_vals_cur);

        n_faces = n_vals_cur - 1;
        val_offset_cur = mr->n_faces_read;

        /* Allocate for first file read */
        if (mr->face_vertices_idx == NULL)
          BFT_MALLOC(mr->face_vertices_idx, n_vals, cs_lnum_t);

        if (val_offset_cur > 0)
          idx_offset_shift = mr->face_vertices_idx[val_offset_cur];

        /* Read data */
        cs_io_set_cs_gnum(&header, pp_in);
        BFT_MALLOC(_g_face_vertices_idx, n_vals_cur+1, cs_gnum_t);
        cs_io_read_index_block(&header, gnum_range_cur[0], gnum_range_cur[1],
                               _g_face_vertices_idx, pp_in);

        /* save start and end values for next read */
        face_vtx_range[1] = _g_face_vertices_idx[n_vals_cur - 1];
        face_vtx_range[0] = _g_face_vertices_idx[0];

        /* Shift cell numbers in case of appended data */
        idx_gnum_shift = _g_face_vertices_idx[0];
        for (ii = 0; ii < n_vals_cur; ii++) {
          cs_gnum_t _face_vtx_idx = _g_face_vertices_idx[ii] - idx_gnum_shift;
          mr->face_vertices_idx[ii + val_offset_cur]
            = _face_vtx_idx + idx_offset_shift;
        }

        BFT_FREE(_g_face_vertices_idx);
      }

      else if (strncmp(header.sec_name, "face_vertices",
                       CS_IO_NAME_LEN) == 0) {

        n_vals = 0;
        n_g_face_connect_size = header.n_vals;

        if (  (cs_gnum_t)(header.n_vals) + mr->n_g_faces_connect_read
            > mr->n_g_face_connect_size)
          bft_error
            (__FILE__, __LINE__, 0,
             _("Section of type <%s> on <%s>\n"
               "has incorrect size (current: %llu, read: %llu, total: %llu."),
             header.sec_name, cs_io_get_name(pp_in),
             (unsigned long long)(header.n_vals),
             (unsigned long long)mr->n_g_faces_connect_read,
             (unsigned long long)mr->n_g_face_connect_size);

        else if (header.n_location_vals != 1)
          bft_error(__FILE__, __LINE__, 0,
                    _("Section of type <%s> on <%s>\n"
                      "has incorrect number of values per location."),
                    header.sec_name, cs_io_get_name(pp_in));

        n_vals_cur = face_vtx_range[1] - face_vtx_range[0];
        val_offset_cur = mr->n_faces_connect_read;

        /* Reallocate for each read, as size of indexed array
           cannot be determined before reading the previous section
           (and is thus not yet known for future files). */
        BFT_REALLOC(mr->face_vertices,
                    mr->n_faces_connect_read + n_vals_cur,
                    cs_gnum_t);

        /* Read data */
        cs_io_set_cs_gnum(&header, pp_in);
        cs_io_read_block(&header, face_vtx_range[0], face_vtx_range[1],
                         mr->face_vertices + val_offset_cur, pp_in);

        /* Shift referenced vertex numbers in case of appended data */
        if (mr->n_g_vertices_read > 0) {
          cs_lnum_t ii;
          for (ii = 0; ii < n_vals_cur; ii++) {
            if (mr->face_vertices[val_offset_cur + ii] != 0)
              mr->face_vertices[val_offset_cur + ii] += mr->n_g_vertices_read;
          }
        }

        mr->n_faces_connect_read += n_vals_cur;
      }

      else if (strncmp(header.sec_name, "vertex_coords",
                       CS_IO_NAME_LEN) == 0) {

        /* Compute range for current file  */
        _data_range(&header,
                    pp_in,
                    mesh->n_g_vertices,
                    mr->n_g_vertices_read,
                    3,
                    0,
                    mr->vertex_bi.gnum_range,
                    gnum_range_cur,
                    &n_g_vertices,
                    &n_vals,
                    &n_vals_cur);

        n_vertices = n_vals_cur / 3;
        val_offset_cur = mr->n_vertices_read * 3;

        /* Allocate for first file read */
        if (mr->vertex_coords == NULL)
          BFT_MALLOC(mr->vertex_coords, n_vals, cs_real_t);

        /* Read data */
        cs_io_assert_cs_real(&header, pp_in);
        cs_io_read_block(&header, gnum_range_cur[0], gnum_range_cur[1],
                         mr->vertex_coords + val_offset_cur, pp_in);

        /* Transform coordinates if necessary */

        if (f->matrix != NULL) {
          cs_gnum_t range_size = gnum_range_cur[1] - gnum_range_cur[0];
          _transform_coords(range_size,
                            mr->vertex_coords + val_offset_cur,
                            f->matrix);
          mesh->modified = 1;
        }
      }

      /* Additional buffers for periodicity */

      else if (strncmp(header.sec_name, "periodicity_type_",
                       strlen("periodicity_type_")) == 0) {

        if (data_read != true || header.n_vals != 1)
          bft_error(__FILE__, __LINE__, 0,
                    _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
        else {
          perio_num = atoi(header.sec_name + strlen("periodicity_type_"));
          n_perio_read = CS_MAX(n_perio_read, perio_num);
          cs_io_read_global(&header, &perio_type, pp_in);
        }

      }
      else if (strncmp(header.sec_name, "periodicity_matrix_",
                       strlen("periodicity_matrix_")) == 0) {

        n_vals = 12; /* 3x4 */
        if (data_read != true || header.n_vals != n_vals)
          bft_error(__FILE__, __LINE__, 0,
                    _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
        else {
          assert(   perio_num
                 == atoi(header.sec_name + strlen("periodicity_matrix_")));
          cs_io_assert_cs_real(&header, pp_in);
          cs_io_read_global(&header, perio_matrix, pp_in);

          /* Add a periodicity to mesh->periodicities */

          if (f->matrix != NULL)
            _transform_perio_matrix(f->matrix, perio_matrix);

          _add_periodicity(mesh,
                           perio_type,
                           perio_num + mr->n_perio_read,
                           perio_matrix);

        }

      }
      else if (strncmp(header.sec_name, "periodicity_faces_",
                       strlen("periodicity_faces_")) == 0) {

        perio_id = atoi(header.sec_name
                        + strlen("periodicity_faces_")) - 1
                        + mr->n_perio_read;

        n_vals = mr->n_g_per_face_couples[perio_id] * 2;

        if (data_read != true || header.n_vals != n_vals)
          bft_error(__FILE__, __LINE__, 0,
                    _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
        else {

          if ((mr->per_face_bi[perio_id]).gnum_range[0] > 0)
            mr->n_per_face_couples[perio_id]
              = (  (mr->per_face_bi[perio_id]).gnum_range[1]
                 - (mr->per_face_bi[perio_id]).gnum_range[0]);
          else
            mr->n_per_face_couples[perio_id] = 0;

          cs_io_set_cs_gnum(&header, pp_in);
          n_vals = mr->n_per_face_couples[perio_id]*2;
          BFT_MALLOC(mr->per_face_couples[perio_id], n_vals, cs_gnum_t);
          assert(header.n_location_vals == 2);
          cs_io_read_block(&header,
                           (mr->per_face_bi[perio_id]).gnum_range[0],
                           (mr->per_face_bi[perio_id]).gnum_range[1],
                           mr->per_face_couples[perio_id],
                           pp_in);

          /* Shift referenced face numbers in case of appended data */
          if (mr->n_g_faces_read > 0) {
            cs_lnum_t ii;
            for (ii = 0; ii < n_vals; ii++) {
              if (mr->per_face_couples[perio_id][ii] != 0)
                mr->per_face_couples[perio_id][ii] += mr->n_g_faces_read;
            }
          }
        }
      }
    }

  } /* End of loop on sections */

  mr->n_perio_read += n_perio_read;
  mr->n_cells_read += n_cells;
  mr->n_faces_read += n_faces;
  mr->n_faces_connect_read += n_face_connect_size;
  mr->n_vertices_read += n_vertices;
  mr->n_g_cells_read += n_g_cells;
  mr->n_g_faces_read += n_g_faces;
  mr->n_g_faces_connect_read += n_g_face_connect_size;
  mr->n_g_vertices_read += n_g_vertices;

  /* Finalize pre-processor input */
  /*------------------------------*/

  f->offset = 0;
  cs_io_finalize(&pp_in);
}

/*----------------------------------------------------------------------------
 * Descend binary tree for the ordering of a mesh's groups.
 *
 * parameters:
 *   mesh   <-- pointer to mesh structure
 *   level  <-- level of the binary tree to descend
 *   n      <-- number of groups in the binary tree to descend
 *   order  <-> ordering array
 *----------------------------------------------------------------------------*/

inline static void
_groups_descend_tree(const cs_mesh_t  *mesh,
                     size_t            level,
                     const size_t      n,
                     int               order[])
{
  size_t i_save, i1, i2, lv_cur;

  i_save = (size_t)(order[level]);

  while (level <= (n/2)) {

    lv_cur = (2*level) + 1;

    if (lv_cur < n - 1) {

      i1 = (size_t)(order[lv_cur+1]);
      i2 = (size_t)(order[lv_cur]);

      if (strcmp(mesh->group_lst + (mesh->group_idx[i1] - 1),
                 mesh->group_lst + (mesh->group_idx[i2] - 1)) > 0)
        lv_cur++;
    }

    if (lv_cur >= n) break;

    i1 = i_save;
    i2 = (size_t)(order[lv_cur]);

    if (strcmp(mesh->group_lst + (mesh->group_idx[i1] - 1),
               mesh->group_lst + (mesh->group_idx[i2] - 1)) >= 0)
      break;

    order[level] = order[lv_cur];
    level = lv_cur;
  }

  order[level] = i_save;
}

/*----------------------------------------------------------------------------
 * Order mesh groups.
 *
 * parameters:
 *   mesh  <-- pointer to mesh structure
 *   order --> pre-allocated ordering table
 *----------------------------------------------------------------------------*/

static void
_order_groups(const cs_mesh_t  *mesh,
              int               order[])
{
  int    o_save;
  size_t i;
  size_t n = mesh->n_groups;

  /* Initialize ordering array */

  for (i = 0; i < n; i++)
    order[i] = i;

  if (n < 2)
    return;

  /* Create binary tree */

  i = (n / 2);
  do {
    i--;
    _groups_descend_tree(mesh, i, n, order);
  } while (i > 0);

  /* Sort binary tree */

  for (i = n - 1; i > 0; i--) {
    o_save   = order[0];
    order[0] = order[i];
    order[i] = o_save;
    _groups_descend_tree(mesh, 0, i, order);
  }
}

/*----------------------------------------------------------------------------
 * Clean mesh group definitions
 *
 * parameters
 *   mesh     <-> pointer to mesh structure
 *----------------------------------------------------------------------------*/

static void
_clean_groups(cs_mesh_t  *mesh)
{
  int i;
  size_t j;
  int n_groups = 0;
  size_t size_tot = 0;
  char *g_prev = NULL, *g_cur = NULL, *g_lst = NULL;
  int *order = NULL, *renum = NULL;

  if (mesh->n_groups < 1)
    return;

  /* Order group names */

  BFT_MALLOC(renum, mesh->n_groups, int);
  BFT_MALLOC(order, mesh->n_groups, int);

  _order_groups(mesh, order);

  /* Build compact group copy */

  BFT_MALLOC(g_lst, mesh->group_idx[mesh->n_groups], char);

  g_cur = mesh->group_lst + (mesh->group_idx[order[0]] - 1);
  g_prev = g_cur;
  n_groups += 1;
  strcpy(g_lst, g_cur);
  size_tot += strlen(g_cur) + 1;
  g_lst[size_tot - 1] = '\0';
  renum[order[0]] = 0;

  for (i = 1; i < mesh->n_groups; i++) {
    g_cur = mesh->group_lst + (mesh->group_idx[order[i]] - 1);
    if (strcmp(g_cur, g_prev) != 0) {
      g_prev = g_cur;
      strcpy(g_lst + size_tot, g_cur);
      n_groups += 1;
      size_tot += strlen(g_cur) + 1;
      g_lst[size_tot - 1] = '\0';
    }
    renum[order[i]] = n_groups - 1;
  }

  BFT_FREE(order);

  BFT_REALLOC(mesh->group_idx, n_groups + 1, cs_int_t);
  BFT_REALLOC(mesh->group_lst, size_tot, char);

  mesh->n_groups = n_groups;
  memcpy(mesh->group_lst, g_lst, size_tot);

  mesh->group_idx[0] = 1;
  for (i = 0; i < mesh->n_groups; i++) {
    j = strlen(mesh->group_lst + (mesh->group_idx[i] - 1)) + 1;
    mesh->group_idx[i + 1] = mesh->group_idx[i] + j;
  }

  BFT_FREE(g_lst);

  /* Now renumber groups in group class description */

  size_tot = mesh->n_families * mesh->n_max_family_items;

  for (j = 0; j < size_tot; j++) {
    int gc_id = mesh->family_item[j];
    if (gc_id < 0)
      mesh->family_item[j] = - renum[-gc_id - 1] - 1;
  }

  BFT_FREE(renum);

  /* Remove empty group if present (it should appear first now) */

  if (mesh->n_groups > 1) {

    if ((mesh->group_idx[1] - mesh->group_idx[0]) == 1) {

      size_t new_lst_size = (  mesh->group_idx[mesh->n_groups]
                             - mesh->group_idx[1]);
      for (i = 0; i < mesh->n_groups; i++)
        mesh->group_idx[i] = mesh->group_idx[i+1] - 1;
      mesh->n_groups -= 1;
      memmove(mesh->group_lst, mesh->group_lst+1, new_lst_size);

      BFT_REALLOC(mesh->group_idx, mesh->n_groups + 1, int);
      BFT_REALLOC(mesh->group_lst, new_lst_size, char);

      for (j = 0; j < size_tot; j++) {
        int gc_id = mesh->family_item[j];
        if (gc_id < 0)
          mesh->family_item[j] += 1;
      }

    }
  }

}

/*============================================================================
 *  Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Query or modification of the option for domain partitioning when no
 * partitioning file is present.
 *
 * This function returns 1 or 2 according to the selected algorithm.
 *
 * Fortran interface :
 *
 * subroutine algdom (iopt)
 * *****************
 *
 * integer          iopt        : <-> : choice of the partitioning base
 *                                        0: query
 *                                        1: initial numbering
 *                                        2: Morton curve (bounding box)
 *                                        3: Morton curve (bounding cube)
 *                                        4: Hilbert curve (bounding box)
 *                                        5: Hilbert curve (bounding cube)
 *----------------------------------------------------------------------------*/

void
CS_PROCF(algdom, ALGDOM)(cs_int_t  *iopt)
{
  *iopt = cs_preprocessor_data_part_choice(*iopt);
}

/*----------------------------------------------------------------------------
 * Read sections from the pre-processor about the dimensions of mesh
 * parameters
 *
 * Fortran Interface:
 *
 * subroutine ledevi(ndim   , nfml  , nprfml, iperio, iperot)
 * *****************
 *
 * integer          ndim        : --> : Spacial dimension (3)
 * integer          nfml        : <-- : Number of families
 * integer          nprfml      : <-- : Number of properties per family
 * integer          iperio      : <-- : Periodicity indicator
 * integer          iperot      : <-- : Number of rotation periodicities
 *----------------------------------------------------------------------------*/

void
CS_PROCF(ledevi, LEDEVI)(const cs_int_t   *ndim,
                         cs_int_t         *nfml,
                         cs_int_t         *nprfml,
                         cs_int_t         *iperio,
                         cs_int_t         *iperot)
{
  int file_id;

  cs_mesh_t  *mesh = cs_glob_mesh;
  _mesh_reader_t *mr = NULL;

  /* Initialize parameter values */

  *nfml = 0;
  *nprfml = 0;

  /* Initialize reading of Preprocessor output */

  _set_default_input_if_needed();

  _cs_glob_mesh_reader = _mesh_reader_create(&_n_mesh_files,
                                             &_mesh_file_info);

  _n_max_mesh_files = 0;

  mr = _cs_glob_mesh_reader;

  for (file_id = 0; file_id < mr->n_files; file_id++)
    _read_dimensions(mesh, mr, file_id);

  /* Return values */

  assert(mesh->dim == *ndim);

  *nfml = mesh->n_families;
  *nprfml = mesh->n_max_family_items;

  if (mesh->n_init_perio > 0)
    *iperio = 1;
  if (mesh->have_rotation_perio > 0)
    *iperot = 1;

  mesh->n_domains = cs_glob_n_ranks;
  mesh->domain_num = cs_glob_rank_id + 1;

  /* Update data in cs_mesh_t structure in serial mode */

  if (cs_glob_n_ranks == 1) {
    mesh->n_cells = mesh->n_g_cells;
    mesh->n_cells_with_ghosts = mesh->n_cells;
    mesh->domain_num = 1;
  }
  else
    mesh->domain_num = cs_glob_rank_id + 1;

  /* Clean group names */

  _clean_groups(mesh);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Query or modification of the option for domain partitioning when no
 * partitioning file is present.
 *
 *  0 : query
 *  1 : based on initial numbering
 *  2 : based on Morton space-filling curve in bounding box
 *  3 : based on Morton space-filling curve in bounding cube
 *  4 : based on Hilbert space-filling curve in bounding box
 *  5 : based on Hilbert space-filling curve in bounding cube (default)
 *
 * choice <-- of partitioning algorithm.
 *
 * returns:
 *   1 to 5 according to the selected algorithm.
 *----------------------------------------------------------------------------*/

int
cs_preprocessor_data_part_choice(int choice)
{
  int retval = 0;

  if (choice < 0 || choice > 5)
    bft_error(__FILE__, __LINE__,0,
              _("The algorithm selection indicator for domain partitioning\n"
                "can take the following values:\n"
                "  1:   partition based on initial numbering\n"
                "  2-5: partition based on space-filling curve\n"
                "and not %d."), choice);

  if (choice == 1)
    _use_sfc = false;
  else if (choice >= 2) {
    _use_sfc = true;
    _sfc_type = choice - 2;
  }

  if (_use_sfc == true)
    retval = _sfc_type + 2;
  else
    retval = 1;

  return retval;
}

/*----------------------------------------------------------------------------
 * Define input mesh file to read.
 *
 * If this function is never called, the default file is read.
 * The first time this function is called,  this default is overriden by the
 * defined file, and all subsequent calls define additional meshes to read.
 *
 * parameters:
 *   file_name       <-- name of file to read
 *   n_group_renames <-- number of groups to rename
 *   group_rename    <-- old (group_rename[i*2]) to new (group_rename[i*2 + 1])
 *                       group names array (size: n_group_renames*2)
 *   transf_matrix   <-- coordinate transformation matrix (or NULL)
 *----------------------------------------------------------------------------*/

void
cs_preprocessor_data_add_file(const char     *file_name,
                              size_t          n_group_renames,
                              const char    **group_rename,
                              const double    transf_matrix[3][4])
{
  size_t  i, l;
  size_t  data_size = 0;
  char  **_old_group_names = NULL, **_new_group_names = NULL;
  _mesh_file_info_t  *f = NULL;

  /* Compute data size */

  data_size = _align_size(strlen(file_name) + 1);

  if (transf_matrix != NULL)
    data_size += _align_size(12*sizeof(double));

  data_size += (_align_size(n_group_renames * sizeof(char *)) * 2);

  for (i = 0; i < n_group_renames; i++) {
    data_size += _align_size(strlen(group_rename[i*2]) + 1);
    if (group_rename[i*2+1] != NULL)
      data_size += _align_size(strlen(group_rename[i*2+1]) + 1);
  }

  /* Allocate data (reallocate mesh file info array f necesary) */

  if (_n_max_mesh_files == 0) {
    _n_max_mesh_files = 1;
    BFT_MALLOC(_mesh_file_info, 1, _mesh_file_info_t);
  }

  if (_n_mesh_files + 1 > _n_max_mesh_files) {
    _n_max_mesh_files *= 2;
    BFT_REALLOC(_mesh_file_info, _n_max_mesh_files, _mesh_file_info_t);
  }

  f = _mesh_file_info + _n_mesh_files;
  _n_mesh_files += 1;

  /* Setup base structeure fields */

  f->offset = 0;
  f->data_size = data_size;
  BFT_MALLOC(f->data, f->data_size, unsigned char);
  memset(f->data, 0, f->data_size);

  /* Setup data */

  data_size = 0;

  l = strlen(file_name) + 1;
  memcpy(f->data, file_name, l);
  f->filename = (const char *)(f->data);

  data_size = _align_size(l);

  if (transf_matrix != NULL) {
    l = 12*sizeof(double);
    memcpy(f->data + data_size, transf_matrix, l);
    f->matrix = (const double *)(f->data + data_size);
    data_size += _align_size(l);
  }
  else
    f->matrix = NULL;

  f->n_group_renames = n_group_renames;
  f->old_group_names = NULL;
  f->new_group_names = NULL;

  if (n_group_renames > 0) {

    _old_group_names = (char **)(f->data + data_size);
    f->old_group_names = (const char *const *)_old_group_names;
    data_size += _align_size(n_group_renames * sizeof(char *));

    _new_group_names = (char **)(f->data + data_size);
    f->new_group_names = (const char *const *)_new_group_names;
    data_size += _align_size(n_group_renames * sizeof(char *));

  }

  for (i = 0; i < n_group_renames; i++) {
    l = strlen(group_rename[i*2]) + 1;
    _old_group_names[i] = (char *)(f->data + data_size);
    memcpy(_old_group_names[i], group_rename[i*2], l);
    data_size += _align_size(l);
    if (group_rename[i*2+1] != NULL) {
      l = strlen(group_rename[i*2+1]) + 1;
      _new_group_names[i] = (char *)(f->data + data_size);
      memcpy(_new_group_names[i], group_rename[i*2+1], l);
      data_size += _align_size(l);
    }
    else
      _new_group_names[i] = NULL;
  }
}

/*----------------------------------------------------------------------------
 * Read pre-processor mesh data and finalize input.
 *
 * At this stage, ghost cells are not generated yet, so the interior
 * face connectivity is not complete near parallel domain or periodic
 * boundaries. Also, isolated faces, if present, are considered to be
 * boundary faces, as they may participate in future mesh joining
 * operations. Their matching cell number will be set to -1.
 * Remaining isolated faces should be removed before completing
 * the mesh structure (for example, just before building ghost cells).
 *
 * parameters:
 *   mesh         <-- pointer to mesh structure
 *   mesh_builder <-- pointer to mesh builder structure
 *----------------------------------------------------------------------------*/

void
cs_preprocessor_data_read_mesh(cs_mesh_t          *mesh,
                               cs_mesh_builder_t  *mesh_builder)
{
  int file_id;

  long  echo = CS_IO_ECHO_OPEN_CLOSE;
  _mesh_reader_t  *mr = _cs_glob_mesh_reader;

  _set_block_ranges(mesh, mr);

  for (file_id = 0; file_id < mr->n_files; file_id++)
    _read_data(file_id, mesh, mr, echo);

  if (mr->n_files > 1)
    mesh->modified = 1;

  /* Read cell rank data if available */

  if (cs_glob_n_ranks > 1 && mr->n_files == 1)
    _read_cell_rank(mesh, mr, echo);

  bft_printf("\n");

  /* Now send data to the correct rank */
  /*-----------------------------------*/

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1)
    _decompose_data_g(mesh,
                      mesh_builder,
                      mr,
                      cs_glob_mpi_comm);

#endif

  if (cs_glob_n_ranks == 1)
    _decompose_data_l(mesh, mesh_builder, mr);

  /* Free temporary memory */

  _mesh_reader_destroy(&mr);
  _cs_glob_mesh_reader = mr;

  /* Remove family duplicates */

  cs_mesh_clean_families(mesh);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
