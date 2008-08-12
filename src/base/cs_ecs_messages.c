/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2008 EDF S.A., France
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
 * Manage the exchange of data between Code_Saturne and the pre-processor
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#if defined(_CS_HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_error.h>
#include <bft_file.h>
#include <bft_mem.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <fvm_periodicity.h>

#include <fvm_block_to_part.h>
#include <fvm_io_num.h>
#include <fvm_interface.h>
#include <fvm_order.h>
#include <fvm_parall.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_mesh.h"
#include "cs_io.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_ecs_messages.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/* Structure used for building mesh structure */
/* ------------------------------------------ */

typedef struct {

  /* Face-related dimensions */

  fvm_gnum_t  n_g_faces;
  fvm_gnum_t  n_g_face_connect_size;

  /* Temporary mesh data */

  int          *cell_rank;

  fvm_gnum_t   *face_cells;
  fvm_lnum_t   *face_vertices_idx;
  fvm_gnum_t   *face_vertices;
  cs_int_t     *cell_gc_id;
  cs_int_t     *face_gc_id;
  cs_real_t    *vertex_coords;

  /* Periodic features */

  int           n_perio;               /* Number of periodicities */
  int          *periodicity_num;       /* Periodicity numbers */
  fvm_lnum_t   *n_per_face_couples;    /* Nb. face couples per periodicity */
  fvm_gnum_t   *n_g_per_face_couples;  /* Global nb. couples per periodicity */

  fvm_gnum_t  **per_face_couples;      /* Periodic face couples list. */

  /* Block ranges for parallel distribution */

  fvm_block_to_part_info_t   cell_bi;     /* Block info for cell data */
  fvm_block_to_part_info_t   face_bi;     /* Block info for face data */
  fvm_block_to_part_info_t   vertex_bi;   /* Block info for vertex data */
  fvm_block_to_part_info_t  *per_face_bi; /* Block info for parallel face
                                             couples */

} _mesh_reader_t;

/*============================================================================
 *  Global variables
 *============================================================================*/

static _mesh_reader_t *_cs_glob_mesh_reader = NULL;

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create an empty mesh reader helper structure.
 *
 * returns:
 *   A pointer to a mesh reader helper structure
 *----------------------------------------------------------------------------*/

static _mesh_reader_t *
_mesh_reader_create(void)
{
  _mesh_reader_t  *mr = NULL;

  BFT_MALLOC(mr, 1, _mesh_reader_t);

  memset(mr, 0, sizeof(_mesh_reader_t));

  mr->n_g_faces = 0;
  mr->n_g_face_connect_size = 0;

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
 * mr <-> pointer to a mesh reader helper
 *
 * returns:
 *  NULL pointer
 *----------------------------------------------------------------------------*/

static void
_mesh_reader_destroy(_mesh_reader_t  **mr)
{
  _mesh_reader_t  *_mr = *mr;

  BFT_FREE(_mr->face_cells);
  BFT_FREE(_mr->face_vertices_idx);
  BFT_FREE(_mr->face_vertices);
  BFT_FREE(_mr->cell_gc_id);
  BFT_FREE(_mr->face_gc_id);
  BFT_FREE(_mr->vertex_coords);

  if (_mr->n_perio > 0) {
    int i;
    for (i = 0; i < _mr->n_perio; i++)
      BFT_FREE(_mr->per_face_couples[i]);
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
 * Parameters:
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
  cs_int_t  i, j, tr_id;
  double  _matrix[3][4];

  fvm_periodicity_type_t _perio_type = perio_type;

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 4; j++)
      _matrix[i][j] = matrix[i][j];
  }

  if (_perio_type == FVM_PERIODICITY_TRANSLATION)
    bft_printf(_(" Ajout de la périodicité %d "
                 "(translation [%10.4e, %10.4e, %10.4e]).\n"),
               (int)perio_num, _matrix[0][3], _matrix[1][3], _matrix[2][3]);

  else if (_perio_type == FVM_PERIODICITY_ROTATION)
    bft_printf(_(" Ajout de la périodicité %d (rotation).\n"),
               (int)perio_num);

  tr_id = fvm_periodicity_add_by_matrix(mesh->periodicity,
                                        perio_num,
                                        FVM_PERIODICITY_ROTATION,
                                        matrix);
}

/*----------------------------------------------------------------------------
 * Set block ranges for parallel reads
 *
 * mesh <-- pointer to mesh structure
 * mr   <-> mesh reader helper
 *----------------------------------------------------------------------------*/

static void
_set_block_ranges(cs_mesh_t       *mesh,
                  _mesh_reader_t  *mr)
{
  int i;

  int rank_id = cs_glob_base_rang;
  int n_ranks = cs_glob_base_nbr;

  /* Always build per_face_range in case of periodicity */

  if (mr->n_perio > 0) {
    BFT_MALLOC(mr->per_face_bi, mr->n_perio, fvm_block_to_part_info_t);
    memset(mr->per_face_bi, 0, sizeof(fvm_block_to_part_info_t)*mr->n_perio);
  }

  /* Set block sizes and ranges (useful for parallel mode) */

  mr->cell_bi = fvm_block_to_part_compute_sizes(rank_id,
                                                n_ranks,
                                                0,
                                                mesh->n_g_cells);

  mr->face_bi = fvm_block_to_part_compute_sizes(rank_id,
                                                n_ranks,
                                                0,
                                                mr->n_g_faces);

  mr->vertex_bi = fvm_block_to_part_compute_sizes(rank_id,
                                                  n_ranks,
                                                  0,
                                                  mesh->n_g_vertices);
  for (i = 0; i < mr->n_perio; i++)
    mr->per_face_bi[i]
      = fvm_block_to_part_compute_sizes(rank_id,
                                        n_ranks,
                                        0,
                                        mr->n_g_per_face_couples[i]);
}

/*----------------------------------------------------------------------------
 * Read cell rank if available
 *
 * mesh <-- pointer to mesh structure
 * mr   <-> mesh reader helper
 * echo <-- echo (verbosity) level
 *----------------------------------------------------------------------------*/

static void
_read_cell_rank(cs_mesh_t       *mesh,
                _mesh_reader_t  *mr,
                size_t           echo)
{
  char file_name[32]; /* more than enough for "domain_number_<n_ranks>" */
  size_t  i;
  cs_io_sec_header_t  header;

  cs_io_t  *rank_pp_in = NULL;
  fvm_lnum_t   n_ranks = 0;
  fvm_gnum_t   n_elts = 0;
  fvm_gnum_t   n_g_cells = 0;

  const char  *unexpected_msg = N_("Message de type <%s> sur <%s>\n"
                                   "inattendu ou de taille incorrecte");

  if (n_ranks == 1)
    return;

#if (_CS_STDC_VERSION < 199901L)
  sprintf(file_name, "domain_number_%d", cs_glob_base_nbr);
#else
  snprintf(file_name, 32, "domain_number_%d", cs_glob_base_nbr);
#endif
  file_name[31] = '\0'; /* Just in case; processor counts would need to be
                           in the exa-range for this to be necessary. */

  /* Test if file exists */

  if (! bft_file_isreg(file_name)) {
    bft_printf(_(" Pas de fichier \"%s\" disponible ;\n"
                 "   on utilisera un découpage de domaines non optimisé.\n"),
               file_name);
    return;
  }

  /* Open file */

  rank_pp_in = cs_io_initialize(file_name,
                                "Domain partitioning, R0",
                                CS_IO_MODE_READ,
                                echo);

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
        cs_io_set_fvm_gnum(&header, rank_pp_in);
        cs_io_read_global(&header, &n_g_cells, rank_pp_in);
        if (n_g_cells != mesh->n_g_cells)
          bft_error(__FILE__, __LINE__, 0,
                    _("Le nombre de cellules indiqué par le fichier\n"
                      "\"%s\" (%lu)\n"
                      "ne correspond pas au nombre à celui du maillage (%lu)."),
                    cs_io_get_name(rank_pp_in),
                    (unsigned long)(n_g_cells),
                    (unsigned long)(mesh->n_g_cells));
      }

    }
    else if (strncmp(header.sec_name, "n_ranks",
                     CS_IO_NAME_LEN) == 0) {

      if (header.n_vals != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name,
                  cs_io_get_name(rank_pp_in));
      else {
        cs_io_set_fvm_lnum(&header, rank_pp_in);
        cs_io_read_global(&header, &n_ranks, rank_pp_in);
        if (n_ranks != cs_glob_base_nbr)
          bft_error(__FILE__, __LINE__, 0,
                    _("Le nombre de rangs indiqué par le fichier\n"
                      "\"%s\" (%d)\n"
                      "ne correspond pas au nombre de rangs courant (%d)."),
                    cs_io_get_name(rank_pp_in), (int)n_ranks,
                    (int)cs_glob_base_nbr);
      }

    }
    else if (strncmp(header.sec_name, "cell:domain number",
                     CS_IO_NAME_LEN) == 0) {

      n_elts = mesh->n_g_cells;
      if (header.n_vals != n_elts)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name,
                  cs_io_get_name(rank_pp_in));
      else {
        cs_io_set_fvm_lnum(&header, rank_pp_in);
        if (mr->cell_bi.gnum_range[0] > 0)
          n_elts = mr->cell_bi.gnum_range[1] - mr->cell_bi.gnum_range[0];
        BFT_MALLOC(mr->cell_rank, n_elts, fvm_lnum_t);
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
  }

  if (rank_pp_in != NULL)
    cs_io_finalize(&rank_pp_in);
}

#if defined(_CS_HAVE_MPI)

/*----------------------------------------------------------------------------
 * Mark faces by type (0 for interior, 1 for exterior faces with outwards
 * pointing normal, 2 for exterior faces with inwards pointing normal,
 * 3 for isolated faces) in parallel mode.
 *
 * The mesh structure is also updated with face counts and connectivity sizes.
 *
 * mesh              <-> pointer to mesh structure
 * n_faces           <-- number of local faces
 * face_ifs          <-- parallel and periodic faces interfaces set
 * face_cell         <-- local face -> cell connectivity
 * face_vertices_idx <-- local face -> vertices index
 * face_type         --> face type marker
 *----------------------------------------------------------------------------*/

static void
_face_type_g(cs_mesh_t                  *mesh,
             fvm_lnum_t                  n_faces,
             const fvm_interface_set_t  *face_ifs,
             const fvm_lnum_t            face_cell[],
             const fvm_lnum_t            face_vertices_idx[],
             char                        face_type[])
{
  fvm_lnum_t i;
  int j;

  const int n_interfaces = fvm_interface_set_size(face_ifs);

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

  for (j = 0; j < n_interfaces; j++) {

    const fvm_interface_t *face_if = fvm_interface_set_get(face_ifs, j);
    fvm_lnum_t face_if_size = fvm_interface_size(face_if);
    const fvm_lnum_t *loc_num = fvm_interface_get_local_num(face_if);

    for (i = 0; i < face_if_size; i++)
      face_type[loc_num[i] - 1] = '\0';

  }

  /* Now count faces of each type */

  mesh->n_i_faces = 0;
  mesh->n_b_faces = 0;
  mesh->i_face_vtx_connect_size = 0;
  mesh->b_face_vtx_connect_size = 0;

  for (i = 0; i < n_faces; i++) {
    fvm_lnum_t n_f_vertices = face_vertices_idx[i+1] - face_vertices_idx[i];
    if (face_type[i] == '\0') {
      mesh->n_i_faces += 1;
      mesh->i_face_vtx_connect_size += n_f_vertices;
    }
    else if (face_type[i] == '\1' || face_type[i] == '\2') {
      mesh->n_b_faces += 1;
      mesh->b_face_vtx_connect_size += n_f_vertices;
    }
  }
}

#endif /* defined(_CS_HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Mark faces by type (0 for interior, 1 for exterior faces with outwards
 * pointing normal, 2 for exterior faces with inwards pointing normal,
 * 3 for isolated faces) in serial mode.
 *
 * The mesh structure is also updated with face counts and connectivity sizes.
 *
 * mesh               <-> pointer to mesh structure
 * n_faces            <-- number of local faces
 * n_periodic_couples <-- number of periodic couples associated with
 *                        each periodic list
 * periodic_couples   <-- array indicating periodic couples (using
 *                        global numberings) for each list
 * face_cell          <-- local face -> cell connectivity
 * face_vertices_idx  <-- local face -> vertices index
 * face_type          --> face type marker
 *----------------------------------------------------------------------------*/

static void
_face_type_l(cs_mesh_t                  *mesh,
             fvm_lnum_t                  n_faces,
             const fvm_lnum_t            n_periodic_couples[],
             const fvm_gnum_t     *const periodic_couples[],
             const fvm_lnum_t            face_cell[],
             const fvm_lnum_t            face_vertices_idx[],
             char                        face_type[])
{
  fvm_lnum_t i;
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

    const fvm_gnum_t *p_couples = periodic_couples[i];

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
    fvm_lnum_t n_f_vertices = face_vertices_idx[i+1] - face_vertices_idx[i];
    if (face_type[i] == '\0') {
      mesh->n_i_faces += 1;
      mesh->i_face_vtx_connect_size += n_f_vertices;
    }
    else if (face_type[i] == '\1' || face_type[i] == '\2') {
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
 * The corresponding arrays in the mesh structure are allocated and
 * defined by this function, and should have been previously empty.
 *
 * mesh      <-> pointer to mesh structure
 * n_faces   <-- number of local faces
 * face_cell <-- local face -> cell connectivity
 * face_type <-- face type marker
 *----------------------------------------------------------------------------*/

static void
_extract_face_cell(cs_mesh_t         *mesh,
                   fvm_lnum_t         n_faces,
                   const fvm_lnum_t   face_cell[],
                   const char         face_type[])
{
  fvm_lnum_t i;

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
  }
}

/*----------------------------------------------------------------------------
 * Build internal and boundary face -> vertices connectivity using a common
 * face -> vertices connectivity and a face type marker.
 *
 * The corresponding arrays in the mesh structure are allocated and
 * defined by this function, and should have been previously empty.
 *
 * mesh              <-> pointer to mesh structure
 * n_faces           <-- number of local faces
 * face_vertices_idx <-- local face -> vertices index
 * face_vertices     <-- local face -> vertices connectivity
 * face_type         <-- face type marker
 *----------------------------------------------------------------------------*/

static void
_extract_face_vertices(cs_mesh_t         *mesh,
                       fvm_lnum_t         n_faces,
                       const fvm_lnum_t   face_vertices_idx[],
                       const fvm_lnum_t   face_vertices[],
                       const char         face_type[])
{
  fvm_lnum_t i;
  size_t j;

  size_t n_i_faces = 0;
  size_t n_b_faces = 0;

  /* Allocate and initialize */

  BFT_MALLOC(mesh->i_face_vtx_idx, mesh->n_i_faces+1, cs_int_t);
  BFT_MALLOC(mesh->i_face_vtx_lst, mesh->i_face_vtx_connect_size, cs_int_t);

  mesh->i_face_vtx_idx[0] = 1;

  if (mesh->n_b_faces > 0) {

    BFT_MALLOC(mesh->b_face_vtx_idx, mesh->n_b_faces+1, cs_int_t);
    BFT_MALLOC(mesh->b_face_vtx_lst, mesh->b_face_vtx_connect_size, cs_int_t);

    mesh->b_face_vtx_idx[0] = 1;

  }

  /* Now copy face -> vertices connectivity */

  for (i = 0; i < n_faces; i++) {

    size_t n_f_vertices = face_vertices_idx[i+1] - face_vertices_idx[i];
    const fvm_lnum_t *_face_vtx = face_vertices + face_vertices_idx[i];

    if (face_type[i] == '\0') {
      fvm_lnum_t *_i_face_vtx =   mesh->i_face_vtx_lst
                                + mesh->i_face_vtx_idx[n_i_faces] - 1;
      for (j = 0; j < n_f_vertices; j++)
        _i_face_vtx[j] = _face_vtx[j];
      mesh->i_face_vtx_idx[n_i_faces + 1] =   mesh->i_face_vtx_idx[n_i_faces]
                                            + n_f_vertices;
      n_i_faces++;
    }

    else if (face_type[i] == '\1') {
      fvm_lnum_t *_b_face_vtx =   mesh->b_face_vtx_lst
                                + mesh->b_face_vtx_idx[n_b_faces] - 1;
      for (j = 0; j < n_f_vertices; j++)
        _b_face_vtx[j] = _face_vtx[j];
      mesh->b_face_vtx_idx[n_b_faces + 1] =   mesh->b_face_vtx_idx[n_b_faces]
                                            + n_f_vertices;
      n_b_faces++;
    }

    else if (face_type[i] == '\2') {
      fvm_lnum_t *_b_face_vtx =   mesh->b_face_vtx_lst
                                + mesh->b_face_vtx_idx[n_b_faces] - 1;
      for (j = 0; j < n_f_vertices; j++)
        _b_face_vtx[j] = _face_vtx[n_f_vertices - j - 1];
      mesh->b_face_vtx_idx[n_b_faces + 1] =   mesh->b_face_vtx_idx[n_b_faces]
                                            + n_f_vertices;
      n_b_faces++;
    }

  }
}

/*----------------------------------------------------------------------------
 * Build internal and boundary face -> global numberings using a common
 * face group class id and a face type marker.
 *
 * The corresponding arrays in the mesh structure are allocated and
 * defined by this function, and should have been previously empty.
 *
 * mesh            <-> pointer to mesh structure
 * n_faces         <-- number of local faces
 * global_face_num <-- global face numbers
 * face_type       <-- face type marker
 *----------------------------------------------------------------------------*/

static void
_extract_face_gnum(cs_mesh_t         *mesh,
                   fvm_lnum_t         n_faces,
                   const fvm_gnum_t   global_face_num[],
                   const char         face_type[])
{
  fvm_lnum_t i;

  size_t n_i_faces = 0;
  size_t n_b_faces = 0;

  fvm_lnum_t *global_i_face = NULL;
  fvm_lnum_t *global_b_face = NULL;

  fvm_io_num_t *tmp_face_num = NULL;

  /* Allocate arrays (including temporary arrays) */

  BFT_MALLOC(mesh->global_i_face_num, mesh->n_i_faces, fvm_gnum_t);
  BFT_MALLOC(mesh->global_b_face_num, mesh->n_b_faces, fvm_gnum_t);

  BFT_MALLOC(global_i_face, mesh->n_i_faces, fvm_lnum_t);
  BFT_MALLOC(global_b_face, mesh->n_b_faces, fvm_lnum_t);

  /* Now build internal and boundary face lists */

  for (i = 0; i < n_faces; i++) {

    if (face_type[i] == '\0')
      global_i_face[n_i_faces++] = i+1;

    else if (face_type[i] == '\1' || face_type[i] == '\2')
      global_b_face[n_b_faces++] = i+1;

  }

  /* Build an I/O numbering on internal faces to compact the global numbering */

  tmp_face_num = fvm_io_num_create(global_i_face,
                                   global_face_num,
                                   n_i_faces,
                                   0);

  memcpy(mesh->global_i_face_num,
         fvm_io_num_get_global_num(tmp_face_num),
         n_i_faces*sizeof(fvm_gnum_t));

  mesh->n_g_i_faces = fvm_io_num_get_global_count(tmp_face_num);

  assert(fvm_io_num_get_local_count(tmp_face_num) == (fvm_lnum_t)n_i_faces);

  tmp_face_num = fvm_io_num_destroy(tmp_face_num);

  /* Build an I/O numbering on boundary faces to compact the global numbering */

  tmp_face_num = fvm_io_num_create(global_b_face,
                                   global_face_num,
                                   n_b_faces,
                                   0);

  if (n_b_faces > 0)
    memcpy(mesh->global_b_face_num,
           fvm_io_num_get_global_num(tmp_face_num),
           n_b_faces*sizeof(fvm_gnum_t));

  mesh->n_g_b_faces = fvm_io_num_get_global_count(tmp_face_num);

  assert(fvm_io_num_get_local_count(tmp_face_num) == (fvm_lnum_t)n_b_faces);

  tmp_face_num = fvm_io_num_destroy(tmp_face_num);

  /* Free remaining temporary arrays */

  BFT_FREE(global_i_face);
  BFT_FREE(global_b_face);
}

/*----------------------------------------------------------------------------
 * Build internal and boundary face -> group class id using a common
 * face group class id and a face type marker.
 *
 * The corresponding arrays in the mesh structure are allocated and
 * defined by this function, and should have been previously empty.
 *
 * mesh       <-> pointer to mesh structure
 * n_faces    <-- number of local faces
 * face_gc_id <-- local face group class id
 * face_type  <-- face type marker
 *----------------------------------------------------------------------------*/

static void
_extract_face_gc_id(cs_mesh_t        *mesh,
                   fvm_lnum_t         n_faces,
                   const fvm_lnum_t   face_gc_id[],
                   const char         face_type[])
{
  fvm_lnum_t i;

  size_t n_i_faces = 0;
  size_t n_b_faces = 0;

  /* Allocate arrays */

  BFT_MALLOC(mesh->i_face_family, mesh->n_i_faces, cs_int_t);
  BFT_MALLOC(mesh->b_face_family, mesh->n_b_faces, cs_int_t);

  /* Now copy face group class (family) id */

  for (i = 0; i < n_faces; i++) {

    if (face_type[i] == '\0')
      mesh->i_face_family[n_i_faces++] = face_gc_id[i];

    else if (face_type[i] == '\1' || face_type[i] == '\2')
      mesh->b_face_family[n_b_faces++] = face_gc_id[i];

  }
}

/*----------------------------------------------------------------------------
 * Re-orient local periodic couples in mesh builder structure.
 * This is probably not necessary, but allows us to build arrays
 * identical to those produced by the preprocessor in version 1.3,
 * so this step may be removed after sufficient testing.
 *
 * mesh_builder      <-> pointer to mesh builder structure
 * n_init_perio      <-- number of initial periodicities
 * i_face_cell       <-- interior face->cell connectivity
 *----------------------------------------------------------------------------*/

static void
_orient_perio_couples(cs_mesh_builder_t  *mb,
                      int                 n_init_perio,
                      const fvm_lnum_t    i_face_cell[])
{
  fvm_lnum_t i;

  const fvm_lnum_t n_couples = mb->per_face_idx[n_init_perio];

  /* In parallel mode */

  if (mb->per_rank_lst != NULL) {

    const int local_rank = cs_glob_base_rang + 1;

    for (i = 0; i < n_couples; i++) {

      if (mb->per_rank_lst[i] == local_rank) {

        fvm_lnum_t inv_sgn = -1;
        fvm_lnum_t face_num_1 = mb->per_face_lst[i*2];
        fvm_lnum_t face_num_2 = mb->per_face_lst[i*2 + 1];
        if (face_num_1 < 0) {
          inv_sgn = 1;
          face_num_1 = -face_num_1;
        }

        if (i_face_cell[(face_num_1-1)*2] == 0) {
          assert(   i_face_cell[(face_num_1-1)*2 + 1] != 0
                 && i_face_cell[(face_num_2-1)*2] != 0
                 && i_face_cell[(face_num_2-1)*2 + 1] == 0);
          mb->per_face_lst[i*2] = face_num_2 * inv_sgn;
          mb->per_face_lst[i*2 + 1] = face_num_1;
        }
      }
    }
  }

  /* In serial mode */

  else { /* if (mb->per_rank_lst == NULL) */

    for (i = 0; i < n_couples; i++) {

      fvm_lnum_t inv_sgn = -1;
      fvm_lnum_t face_num_1 = mb->per_face_lst[i*2];
      fvm_lnum_t face_num_2 = mb->per_face_lst[i*2 + 1];
      if (face_num_1 < 0) {
        inv_sgn = 1;
        face_num_1 = -face_num_1;
      }

      if (i_face_cell[(face_num_1-1)*2] == 0) {
        mb->per_face_lst[i*2] = face_num_2 * inv_sgn;
        mb->per_face_lst[i*2 + 1] = face_num_1;
      }
    }
  }
}

#if defined(_CS_HAVE_MPI)

/*----------------------------------------------------------------------------
 * Extract periodic face connectivity information for mesh builder when
 * running in parallel mode.
 *
 * mesh_builder      <-> pointer to mesh builder structure
 * n_init_perio      <-- number of initial periodicities
 * n_faces           <-- number of local faces
 * face_ifs          <-- parallel and periodic faces interfaces set
 * face_type         <-- face type marker
 *----------------------------------------------------------------------------*/

static void
_extract_periodic_faces_g(cs_mesh_builder_t          *mb,
                          int                         n_init_perio,
                          fvm_lnum_t                  n_faces,
                          const fvm_interface_set_t  *face_ifs,
                          const char                  face_type[])
{
  fvm_lnum_t i;
  int j;

  fvm_lnum_t   i_face_count = 0;
  fvm_lnum_t  *i_face_id = NULL;
  fvm_lnum_t  *per_face_count = NULL;
  fvm_lnum_t  *if_index = NULL;
  fvm_lnum_t  *send_num = NULL, *recv_num = NULL;

  const int n_interfaces = fvm_interface_set_size(face_ifs);
  const fvm_lnum_t tr_index_size = n_init_perio*2 + 2;

  /* Allocate arrays in mesh builder (initializing per_face_idx) */

  BFT_MALLOC(mb->per_face_idx, n_init_perio + 1, cs_int_t);

  for (i = 0; i < n_init_perio + 1; i++)
    mb->per_face_idx[i] = 0;

  for (j = 0; j < n_interfaces; j++) {

    const fvm_interface_t *face_if = fvm_interface_set_get(face_ifs, j);
    const fvm_lnum_t *tr_index = fvm_interface_get_tr_index(face_if);
    const int distant_rank = fvm_interface_rank(face_if);

    assert(fvm_interface_get_tr_index_size(face_if) == tr_index_size);

    /* Only count 1 transformation direction when corresponding
       faces are on the same rank (in which case they appear
       once per oposing direction transforms) */

    for (i = 1; i < tr_index_size-1; i++) {
      if ((distant_rank != cs_glob_base_rang) || (i%2 == 1))
        mb->per_face_idx[(i-1)/2 + 1] += tr_index[i+1] - tr_index[i];
    }
  }

  mb->per_face_idx[0] = 0;
  for (i = 1; i < n_init_perio+1; i++)
    mb->per_face_idx[i] += mb->per_face_idx[i-1];

  BFT_MALLOC(mb->per_face_lst, mb->per_face_idx[n_init_perio] * 2, cs_int_t);
  BFT_MALLOC(mb->per_rank_lst, mb->per_face_idx[n_init_perio], cs_int_t);

  /* Build face renumbering */

  BFT_MALLOC(i_face_id, n_faces, fvm_lnum_t);

  for (i = 0; i < n_faces; i++) {
    if (face_type[i] == '\0')
      i_face_id[i] = i_face_count++;
    else
      i_face_id[i] = -1;
  }

  /* Copy periodic interface arrays and renumber them */

  BFT_MALLOC(if_index, n_interfaces + 1, fvm_lnum_t);
  if_index[0] = 0;

  for (j = 0; j < n_interfaces; j++) {
    const fvm_interface_t *face_if = fvm_interface_set_get(face_ifs, j);
    const fvm_lnum_t *tr_index = fvm_interface_get_tr_index(face_if);
    if_index[j+1] = if_index[j] + tr_index[tr_index_size - 1] - tr_index[1];
  }

  BFT_MALLOC(send_num, if_index[n_interfaces], fvm_lnum_t);
  BFT_MALLOC(recv_num, if_index[n_interfaces], fvm_lnum_t);

  for (j = 0; j < n_interfaces; j++) {

    fvm_lnum_t k, l;

    const fvm_lnum_t start_id = if_index[j];
    const fvm_lnum_t end_id = if_index[j+1];
    const fvm_interface_t *face_if = fvm_interface_set_get(face_ifs, j);
    const fvm_lnum_t *tr_index = fvm_interface_get_tr_index(face_if);
    const fvm_lnum_t *loc_num = fvm_interface_get_local_num(face_if);
    const int distant_rank = fvm_interface_rank(face_if);

    for (k = start_id, l = tr_index[1]; k < end_id; k++, l++)
      send_num[k] = i_face_id[loc_num[l] - 1] + 1;

    if (distant_rank == cs_glob_base_rang) {
      const fvm_lnum_t *dist_num = fvm_interface_get_distant_num(face_if);
      for (k = start_id, l = tr_index[1]; k < end_id; k++, l++)
        recv_num[k] = i_face_id[dist_num[l] - 1] + 1;
    }
  }

  BFT_FREE(i_face_id);

  /* Exchange local face numbers */

  {
    MPI_Request  *request = NULL;
    MPI_Status  *status  = NULL;

    int request_count = 0;

    BFT_MALLOC(request, n_interfaces*2, MPI_Request);
    BFT_MALLOC(status, n_interfaces*2, MPI_Status);

    for (j = 0; j < n_interfaces; j++) {
      const fvm_interface_t *face_if = fvm_interface_set_get(face_ifs, j);
      int distant_rank = fvm_interface_rank(face_if);
      if (distant_rank != cs_glob_base_rang)
        MPI_Irecv(recv_num + if_index[j],
                  if_index[j+1] - if_index[j],
                  FVM_MPI_LNUM,
                  distant_rank,
                  distant_rank,
                  cs_glob_base_mpi_comm,
                  &(request[request_count++]));
    }

    for (j = 0; j < n_interfaces; j++) {
      const fvm_interface_t *face_if = fvm_interface_set_get(face_ifs, j);
      int distant_rank = fvm_interface_rank(face_if);
      if (distant_rank != cs_glob_base_rang)
        MPI_Isend(send_num + if_index[j],
                  if_index[j+1] - if_index[j],
                  FVM_MPI_LNUM,
                  distant_rank,
                  (int)cs_glob_base_rang,
                  cs_glob_base_mpi_comm,
                  &(request[request_count++]));
    }

    MPI_Waitall(request_count, request, status);

    BFT_FREE(request);
    BFT_FREE(status);
  }

  /* Copy new interface information to mesh builder */

  BFT_MALLOC(per_face_count, n_init_perio, fvm_lnum_t);
  for (i = 0; i < n_init_perio; i++)
    per_face_count[i] = 0;

  for (j = 0; j < n_interfaces; j++) {

    const fvm_interface_t *face_if = fvm_interface_set_get(face_ifs, j);
    const int distant_rank = fvm_interface_rank(face_if);
    const fvm_lnum_t *tr_index = fvm_interface_get_tr_index(face_if);

    for (i = 1; i < tr_index_size - 1; i++) {

      if ((distant_rank != cs_glob_base_rang) || (i%2 == 1)) {

        fvm_lnum_t k, l;

        int perio_id = (i-1)/2;
        int perio_sgn = (i%2)*2 - 1; /* 1 for odd, -1 for even */
        fvm_lnum_t tr_start_id = if_index[j] + (tr_index[i] - tr_index[1]);
        fvm_lnum_t tr_end_id = if_index[j] + (tr_index[i+1] - tr_index[1]);

        for (k = tr_start_id; k < tr_end_id; k++) {
          l = mb->per_face_idx[perio_id] + per_face_count[perio_id];
          mb->per_face_lst[l*2]     = send_num[k]*perio_sgn;
          mb->per_face_lst[l*2 + 1] = recv_num[k];
          mb->per_rank_lst[l] = distant_rank + 1;
          per_face_count[perio_id] += 1;
        }
      }
    }
  }

  BFT_FREE(per_face_count);
  BFT_FREE(recv_num);
  BFT_FREE(send_num);
  BFT_FREE(if_index);
}

#endif /* defined(_CS_HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Extract periodic face connectivity information for mesh builder when
 * running in serial mode.
 *
 * mesh_builder       <-> pointer to mesh builder structure
 * n_init_perio       <-- number of initial periodicities
 * n_faces            <-- number of local faces
 * n_periodic_couples <-- number of periodic couples associated with
 *                        each periodic list
 * periodic_couples   <-- array indicating periodic couples (using
 *                        global numberings) for each list
 * face_type          <-- face type marker
 *----------------------------------------------------------------------------*/

static void
_extract_periodic_faces_l(cs_mesh_builder_t        *mb,
                          int                       n_init_perio,
                          fvm_lnum_t                n_faces,
                          const fvm_lnum_t          n_periodic_couples[],
                          const fvm_gnum_t   *const periodic_couples[],
                          const char                face_type[])
{
  int i;
  fvm_lnum_t j;

  fvm_lnum_t   i_face_count = 0;
  fvm_lnum_t  *i_face_id = NULL;

  /* Allocate arrays in mesh builder (initializing per_face_idx) */

  BFT_MALLOC(mb->per_face_idx, n_init_perio + 1, cs_int_t);

  mb->per_face_idx[0] = 0;
  for (i = 0; i < n_init_perio; i++)
    mb->per_face_idx[i+1] = mb->per_face_idx[i] + n_periodic_couples[i];

  BFT_MALLOC(mb->per_face_lst, mb->per_face_idx[n_init_perio] * 2, cs_int_t);

  /* Build face renumbering */

  BFT_MALLOC(i_face_id, n_faces, fvm_lnum_t);

  for (i = 0; i < n_faces; i++) {
    if (face_type[i] == '\0')
      i_face_id[i] = i_face_count++;
    else
      i_face_id[i] = -1;
  }

  /* Copy new interface information to mesh builder */

  for (i = 0; i < n_init_perio; i++) {

    const fvm_gnum_t *p_couples = periodic_couples[i];

    for (j = 0; j < n_periodic_couples[i]; j++) {

      fvm_lnum_t k = mb->per_face_idx[i] + j;

      mb->per_face_lst[k*2]   = i_face_id[p_couples[j*2] - 1] + 1;
      mb->per_face_lst[k*2+1] = i_face_id[p_couples[j*2+1] - 1] + 1;

    }

  }

  BFT_FREE(i_face_id);
}

#if defined(_CS_HAVE_MPI)

/*----------------------------------------------------------------------------
 * Organize data read by blocks in parallel and build most mesh structures.
 *
 * mesh         <-- pointer to mesh structure
 * mesh_builder <-- pointer to mesh builder structure
 * mr           <-> pointer to mesh reader helper structure
 *----------------------------------------------------------------------------*/

static void
_decompose_data_g(cs_mesh_t          *mesh,
                  cs_mesh_builder_t  *mesh_builder,
                  _mesh_reader_t     *mr,
                  MPI_Comm            comm)
{
  fvm_lnum_t i;
  int n_ranks = 0;

  fvm_datatype_t lnum_type = (sizeof(fvm_lnum_t) == 8) ? FVM_INT64 : FVM_INT32;
  fvm_datatype_t gnum_type = (sizeof(fvm_gnum_t) == 8) ? FVM_UINT64 : FVM_UINT32;
  fvm_datatype_t real_type = (sizeof(cs_real_t) == 8) ? FVM_DOUBLE : FVM_FLOAT;

  fvm_lnum_t _n_faces = 0;
  fvm_gnum_t cell_block_size = 0;
  fvm_gnum_t face_block_size = 0;
  fvm_gnum_t vertex_block_size = 0;

  fvm_gnum_t *_face_num = NULL;
  fvm_gnum_t *_face_gcells = NULL;
  fvm_gnum_t *_face_gvertices = NULL;

  fvm_lnum_t *_face_cells = NULL;
  fvm_lnum_t *_face_gc_id = NULL;
  fvm_lnum_t *_face_vertices_idx = NULL;
  fvm_lnum_t *_face_vertices = NULL;

  char *face_type = NULL;
  fvm_interface_set_t *face_ifs = NULL;

  fvm_block_to_part_t *d = NULL;

  /* Initialization */

  MPI_Comm_size(comm, &n_ranks);

  cell_block_size = mesh->n_g_cells / n_ranks;
  face_block_size = mr->n_g_faces / n_ranks;
  vertex_block_size = mesh->n_g_vertices / n_ranks;

  assert((sizeof(fvm_lnum_t) == 4) || (sizeof(fvm_lnum_t) == 8));

  /* Different handling of cells depending on whether decomposition
     data is available or not. */

  if (mr->cell_rank != NULL) {

    d = fvm_block_to_part_create_by_rank(comm,
                                         mr->cell_bi,
                                         mr->cell_rank);

    mesh->n_cells = fvm_block_to_part_get_n_part_ents(d);

    BFT_MALLOC(mesh->cell_family, mesh->n_cells, fvm_lnum_t);

    fvm_block_to_part_copy_array(d,
                                 lnum_type,
                                 1,
                                 mr->cell_gc_id,
                                 mesh->cell_family);

    BFT_FREE(mr->cell_gc_id);

    mesh->global_cell_num = fvm_block_to_part_transfer_gnum(d);

    fvm_block_to_part_destroy(&d);

  }
  else {

    mesh->n_cells = mr->cell_bi.gnum_range[1] - mr->cell_bi.gnum_range[0];

    BFT_MALLOC(mesh->global_cell_num, mesh->n_cells, fvm_gnum_t);

    for (i = 0; i < mesh->n_cells; i++)
      mesh->global_cell_num[i] = mr->cell_bi.gnum_range[0] + i;

    mesh->cell_family = mr->cell_gc_id;
    mr->cell_gc_id = NULL;
  }

  /* Distribute faces */
  /*------------------*/

  d = fvm_block_to_part_create_by_adj_s(comm,
                                        mr->face_bi,
                                        mr->cell_bi,
                                        2,
                                        mr->face_cells,
                                        mr->cell_rank);

  BFT_FREE(mr->cell_rank); /* Not needed anymore */

  _n_faces = fvm_block_to_part_get_n_part_ents(d);

  BFT_MALLOC(_face_gcells, _n_faces*2, fvm_gnum_t);

  /* Face -> cell connectivity */

  fvm_block_to_part_copy_array(d,
                               gnum_type,
                               2,
                               mr->face_cells,
                               _face_gcells);

  BFT_FREE(mr->face_cells);

  /* Now convert face -> cell connectivity to local cell numbers */

  BFT_MALLOC(_face_cells, _n_faces*2, fvm_lnum_t);

  fvm_block_to_part_global_to_local(_n_faces*2,
                                    1,
                                    mesh->n_cells,
                                    mesh->global_cell_num,
                                    _face_gcells,
                                    _face_cells);

  BFT_FREE(_face_gcells);

  /* Face family */

  BFT_MALLOC(_face_gc_id, _n_faces, fvm_lnum_t);

  fvm_block_to_part_copy_array(d,
                               lnum_type,
                               1,
                               mr->face_gc_id,
                               _face_gc_id);

  BFT_FREE(mr->face_gc_id);

  /* Face connectivity */

  BFT_MALLOC(_face_vertices_idx, _n_faces + 1, fvm_lnum_t);

  fvm_block_to_part_copy_index(d,
                               mr->face_vertices_idx,
                               _face_vertices_idx);

  BFT_MALLOC(_face_gvertices, _face_vertices_idx[_n_faces], fvm_gnum_t);

  fvm_block_to_part_copy_indexed(d,
                                 gnum_type,
                                 mr->face_vertices_idx,
                                 mr->face_vertices,
                                 _face_vertices_idx,
                                 _face_gvertices);

  BFT_FREE(mr->face_vertices_idx);
  BFT_FREE(mr->face_vertices);

  _face_num = fvm_block_to_part_transfer_gnum(d);

  fvm_block_to_part_destroy(&d);

  /* Vertices */

  d = fvm_block_to_part_create_adj(comm,
                                   mr->vertex_bi,
                                   _face_vertices_idx[_n_faces],
                                   _face_gvertices);

  mesh->n_vertices = fvm_block_to_part_get_n_part_ents(d);

  BFT_MALLOC(mesh->vtx_coord, mesh->n_vertices*3, cs_real_t);

  fvm_block_to_part_copy_array(d,
                               real_type,
                               3,
                               mr->vertex_coords,
                               mesh->vtx_coord);

  BFT_FREE(mr->vertex_coords);

  mesh->global_vtx_num = fvm_block_to_part_transfer_gnum(d);

  fvm_block_to_part_destroy(&d);

  /* Now convert face -> vertex connectivity to local vertex numbers */

  BFT_MALLOC(_face_vertices, _face_vertices_idx[_n_faces], fvm_lnum_t);

  fvm_block_to_part_global_to_local(_face_vertices_idx[_n_faces],
                                    1,
                                    mesh->n_vertices,
                                    mesh->global_vtx_num,
                                    _face_gvertices,
                                    _face_vertices);

  BFT_FREE(_face_gvertices);

  /* In case of periodicity, build an fvm_interface so as to obtain
     periodic face correspondants in local numbering (periodic couples
     need not be defined by the ranks owning one of the 2 members
     for the interface to be built correctly). */

  face_ifs
    = fvm_interface_set_create(_n_faces,
                               NULL,
                               _face_num,
                               mesh->periodicity,
                               mr->n_perio,
                               mr->periodicity_num,
                               mr->n_per_face_couples,
                               (const fvm_gnum_t **const)mr->per_face_couples);

  /* We may now separate interior from boundary faces */

  BFT_MALLOC(face_type, _n_faces, char);

  _face_type_g(mesh,
               _n_faces,
               face_ifs,
               _face_cells,
               _face_vertices_idx,
               face_type);

  _extract_face_cell(mesh, _n_faces, _face_cells, face_type);

  BFT_FREE(_face_cells);

  if (mr->n_perio > 0) {

    _extract_periodic_faces_g(mesh_builder,
                              mesh->n_init_perio,
                              _n_faces,
                              face_ifs,
                              face_type);

    _orient_perio_couples(mesh_builder,
                          mesh->n_init_perio,
                          mesh->i_face_cells);

  }

  face_ifs = fvm_interface_set_destroy(face_ifs);

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

  _extract_face_gc_id(mesh,
                      _n_faces,
                      _face_gc_id,
                      face_type);

  BFT_FREE(_face_gc_id);

  BFT_FREE(face_type);
}

#endif /* defined(_CS_HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Organize data read locally and build most mesh structures
 *
 * mesh         <-- pointer to mesh structure
 * mesh_builder <-- pointer to mesh builder structure
 * mr           <-> pointer to mesh reader helper structure
 *----------------------------------------------------------------------------*/

static void
_decompose_data_l(cs_mesh_t          *mesh,
                  cs_mesh_builder_t  *mesh_builder,
                  _mesh_reader_t     *mr)
{
  fvm_lnum_t i;

  fvm_lnum_t _n_faces = 0;

  fvm_lnum_t *_face_cells = NULL;
  fvm_lnum_t *_face_vertices_idx = NULL;
  fvm_lnum_t *_face_vertices = NULL;

  char *face_type = NULL;

  /* Initialization */

  assert((sizeof(fvm_lnum_t) == 4) || (sizeof(fvm_lnum_t) == 8));

  mesh->n_cells = mr->cell_bi.gnum_range[1] - 1;

  /* Cell families are already of the correct type,
     so they can simply be moved */

  mesh->cell_family = mr->cell_gc_id;
  mr->cell_gc_id = NULL;

  /* Build faces */
  /*-------------*/

  _n_faces = mr->face_bi.gnum_range[1] - 1;

  /* Now copy face -> cell connectivity to local cell numbers */

  BFT_MALLOC(_face_cells, _n_faces*2, fvm_lnum_t);

  for (i = 0; i < _n_faces; i++) {
    _face_cells[i*2] = mr->face_cells[i*2];
    _face_cells[i*2 + 1] = mr->face_cells[i*2 + 1];
  }

  BFT_FREE(mr->face_cells);

  /* Face connectivity */

  BFT_MALLOC(_face_vertices_idx, _n_faces + 1, fvm_lnum_t);

  for (i = 0; i < _n_faces+1; i++)
    _face_vertices_idx[i] = mr->face_vertices_idx[i];

  BFT_FREE(mr->face_vertices_idx);

  BFT_MALLOC(_face_vertices, _face_vertices_idx[_n_faces], fvm_lnum_t);

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
               (const fvm_gnum_t **const)mr->per_face_couples,
               _face_cells,
               _face_vertices_idx,
               face_type);

  _extract_face_cell(mesh, _n_faces, _face_cells, face_type);

  BFT_FREE(_face_cells);

  if (mr->n_perio > 0) {

    _extract_periodic_faces_l(mesh_builder,
                              mesh->n_init_perio,
                              _n_faces,
                              mr->n_per_face_couples,
                              (const fvm_gnum_t **const)mr->per_face_couples,
                              face_type);

    _orient_perio_couples(mesh_builder,
                          mesh->n_init_perio,
                          mesh->i_face_cells);

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

/*============================================================================
 *  Public functions definition for API Fortran
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Receive messages from the pre-processor about the dimensions of mesh
 * parameters
 *
 * FORTRAN Interface:
 *
 * SUBROUTINE LEDEVI(NOMRUB, TYPENT, NBRENT, TABENT)
 * *****************
 *
 * INTEGER          NDIM        : <-- : Dimension de l'espace (3)
 * INTEGER          NFML        : <-- : Nombre de familles des faces de bord
 * INTEGER          NPRFML      : <-- : Nombre de propriétés max par famille
 * INTEGER          IPERIO      : <-- : Indicateur de périodicité
 * INTEGER          IPEROT      : <-- : Nombre de périodicités de rotation
 *----------------------------------------------------------------------------*/

void CS_PROCF(ledevi, LEDEVI)
(
 cs_int_t   *const ndim,    /* <-- dimension de l'espace                      */
 cs_int_t   *const nfml,    /* <-- nombre de familles des faces de bord       */
 cs_int_t   *const nprfml,  /* <-- nombre de propriétés max par famille       */
 cs_int_t   *const iperio,  /* <-- indicateur de périodicité                  */
 cs_int_t   *const iperot   /* <-- nombre de périodicités de rotation         */
)
{
  cs_int_t  i;
  cs_io_sec_header_t  header;

  fvm_gnum_t n_elts = 0;
  cs_bool_t  dim_read = false;
  cs_bool_t  end_read = false;
  cs_io_t  *pp_in = cs_glob_pp_io;
  cs_mesh_t  *mesh = cs_glob_mesh;
  _mesh_reader_t *mr = NULL;

  const char  *unexpected_msg = N_("Message de type <%s> sur <%s>\n"
                                  "inattendu ou de taille incorrecte");

  /* Initialize parameter values */

  *ndim = 3;
  *nfml = 0;
  *nprfml = 0;

  mr = _mesh_reader_create();

  _cs_glob_mesh_reader = mr;

  /* Loop on read sections */

  while (end_read == false) {

    /* Receive headers and clen header names */

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

    else if (strncmp(header.sec_name, "ndim",
                     CS_IO_NAME_LEN) == 0) {

      if (dim_read != true || header.n_vals != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else
        cs_io_read_global(&header, (void *) &(mesh->dim), pp_in);

    }
    else if (strncmp(header.sec_name, "n_cells",
                     CS_IO_NAME_LEN) == 0) {

      if (dim_read != true || header.n_vals != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else {
        cs_io_set_fvm_gnum(&header, pp_in);
        cs_io_read_global(&header, &(mesh->n_g_cells), pp_in);
      }

    }
    else if (strncmp(header.sec_name, "n_faces",
                     CS_IO_NAME_LEN) == 0) {

      if (dim_read != true || header.n_vals != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else {
        cs_io_set_fvm_gnum(&header, pp_in);
        cs_io_read_global(&header, &(mr->n_g_faces), pp_in);
      }

    }
    else if (strncmp(header.sec_name, "n_vertices",
                     CS_IO_NAME_LEN) == 0) {

      if (dim_read != true || header.n_vals != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else {
        cs_io_set_fvm_gnum(&header, pp_in);
        cs_io_read_global(&header, &(mesh->n_g_vertices), pp_in);
      }

    }
    else if (strncmp(header.sec_name, "face_vertices_size",
                     CS_IO_NAME_LEN) == 0) {

      if (dim_read != true || header.n_vals != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else {
        cs_io_set_fvm_gnum(&header, pp_in);
        cs_io_read_global(&header, &(mr->n_g_face_connect_size), pp_in);
      }

    }
    else if (strncmp(header.sec_name, "n_group_classes",
                     CS_IO_NAME_LEN) == 0) {
      if (dim_read != true || header.n_vals != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else
        cs_io_read_global(&header, (void *) &(mesh->n_families), pp_in);

    }
    else if (strncmp(header.sec_name, "n_group_class_props_max",
                     CS_IO_NAME_LEN) == 0) {

      if (dim_read != true || header.n_vals != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else
        cs_io_read_global(&header,
                          (void *) &(mesh->n_max_family_items), pp_in);

    }
    else if (strncmp(header.sec_name, "n_groups",
                     CS_IO_NAME_LEN) == 0) {

      if (dim_read != true || header.n_vals != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else
        cs_io_read_global(&header, (void *) &(mesh->n_groups), pp_in);

    }
    else if (strncmp(header.sec_name, "group_name_index",
                     CS_IO_NAME_LEN) == 0) {

      if ((cs_int_t)header.n_vals != mesh->n_groups + 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else {
        BFT_MALLOC(mesh->group_idx, mesh->n_groups + 1, cs_int_t);
        cs_io_read_global(&header, (void *) mesh->group_idx, pp_in);
      }

    }
    else if (strncmp(header.sec_name, "group_name",
                     CS_IO_NAME_LEN) == 0) {

      if (   mesh->group_idx == NULL
          || (cs_int_t)header.n_vals != mesh->group_idx[mesh->n_groups] - 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else {
        BFT_MALLOC(mesh->group_lst, header.n_vals, char);
        cs_io_read_global(&header, (void *) mesh->group_lst, pp_in);
      }

    }
    else if (   strncmp(header.sec_name, "group_class_properties",
                        CS_IO_NAME_LEN) == 0
             || strncmp(header.sec_name, "iprfml",
                        CS_IO_NAME_LEN) == 0) {

      n_elts = mesh->n_families * mesh->n_max_family_items;
      if (dim_read != true || header.n_vals != n_elts)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else {
        BFT_MALLOC(mesh->family_item, n_elts, cs_int_t);
        cs_io_read_global(&header, (void *) mesh->family_item, pp_in);
      }

    }

    /* Additional messages for periodicity. Dimensions for periodic ghost
       cells have been received before. Here we allocate parameter list
       for periodicities and coupled face list for halo builder. */

    else if (strncmp(header.sec_name, "n_periodic_directions",
                     CS_IO_NAME_LEN) == 0) {

      if (dim_read != true || header.n_vals != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else {
        cs_io_read_global(&header, (void *) &(mesh->n_init_perio), pp_in);

        assert(mesh->n_init_perio > 0);

        *iperio = 1;
        mesh->periodicity = fvm_periodicity_create(0.001);

        BFT_MALLOC(mr->periodicity_num, mesh->n_init_perio, int);
        BFT_MALLOC(mr->n_per_face_couples, mesh->n_init_perio, fvm_lnum_t);
        BFT_MALLOC(mr->n_g_per_face_couples, mesh->n_init_perio, fvm_gnum_t);
        BFT_MALLOC(mr->per_face_couples, mesh->n_init_perio, fvm_gnum_t *);

        mr->n_perio = mesh->n_init_perio;

        for (i = 0; i < mesh->n_init_perio; i++) {
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
        cs_io_read_global(&header, (void *) iperot, pp_in);
        if (*iperot > 0)
          mesh->have_rotation_perio = 1;
      }

    }
    else if (strncmp(header.sec_name, "n_periodic_faces",
                     CS_IO_NAME_LEN) == 0) {

      if ((cs_int_t)header.n_vals != mesh->n_init_perio)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else {
        cs_io_set_fvm_gnum(&header, pp_in);
        cs_io_read_global(&header, mr->n_g_per_face_couples, pp_in);
        for (i = 0; i < mesh->n_init_perio; i++)
          mr->n_g_per_face_couples[i] /= 2;
      }

    }
    else
      bft_error(__FILE__, __LINE__, 0,
                _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));

  } /* End of test on headers */

  /* Return values */

  *ndim = mesh->dim;
  *nfml = mesh->n_families;
  *nprfml = mesh->n_max_family_items;

  mesh->n_domains = cs_glob_base_nbr;
  mesh->domain_num = cs_glob_base_rang + 1;

  /* Update data in cs_mesh_t structure in serial mode */

  if (cs_glob_base_nbr == 1) {
    mesh->n_cells = mesh->n_g_cells;
    mesh->n_cells_with_ghosts = mesh->n_cells;
    mesh->domain_num = 1;
  }
  else
    mesh->domain_num = cs_glob_base_rang + 1;
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Read data from the pre-processor and finalize pre-processor input.
 *
 * mesh         <-- pointer to mesh structure
 * mesh_builder <-- pointer to mesh builder structure
 *
 * returns:
 *----------------------------------------------------------------------------*/

void
cs_ecs_messages_read_data(cs_mesh_t          *mesh,
                          cs_mesh_builder_t  *mesh_builder)
{
  cs_int_t  perio_id, perio_type;
  cs_io_sec_header_t  header;

  cs_real_t  perio_matrix[3][4];

  cs_int_t  perio_num = -1;
  fvm_gnum_t n_elts = 0;
  fvm_gnum_t face_vertices_idx_shift = 0;
  cs_bool_t  end_read = false;
  cs_bool_t  data_read = false;
  size_t  echo = 0;
  cs_io_t  *pp_in = cs_glob_pp_io;
  _mesh_reader_t  *mr = _cs_glob_mesh_reader;

  const char  *unexpected_msg = N_("Section de type <%s> sur <%s>\n"
                                   "inattendue ou de taille incorrecte");

  echo = cs_io_get_echo(pp_in);

  _set_block_ranges(mesh, mr);

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

    else if (strncmp(header.sec_name, "face_cells",
                     CS_IO_NAME_LEN) == 0) {

      n_elts = mr->n_g_faces * 2;
      if (data_read != true || header.n_vals != n_elts)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else {
        cs_io_set_fvm_gnum(&header, pp_in);
        if (mr->face_bi.gnum_range[0] > 0)
          n_elts = (mr->face_bi.gnum_range[1] - mr->face_bi.gnum_range[0])*2;
        BFT_MALLOC(mr->face_cells, n_elts, fvm_gnum_t);
        cs_io_read_block(&header,
                         mr->face_bi.gnum_range[0]*2 -1,
                         mr->face_bi.gnum_range[1]*2 -1,
                         mr->face_cells, pp_in);
      }

    }
    else if (strncmp(header.sec_name, "cell_group_class_id",
                     CS_IO_NAME_LEN) == 0) {

      n_elts = mesh->n_g_cells;
      if (data_read != true || header.n_vals != n_elts)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else {
        cs_io_set_fvm_lnum(&header, pp_in);
        if (mr->cell_bi.gnum_range[0] > 0)
          n_elts = mr->cell_bi.gnum_range[1] - mr->cell_bi.gnum_range[0];
        BFT_MALLOC(mr->cell_gc_id, n_elts, cs_int_t);
        cs_io_read_block(&header,
                         mr->cell_bi.gnum_range[0],
                         mr->cell_bi.gnum_range[1],
                         mr->cell_gc_id, pp_in);
      }

    }
    else if (strncmp(header.sec_name, "face_group_class_id",
                     CS_IO_NAME_LEN) == 0) {

      n_elts = mr->n_g_faces;
      if (data_read != true || header.n_vals != n_elts)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else {
        cs_io_set_fvm_lnum(&header, pp_in);
        if (mr->face_bi.gnum_range[0] > 0)
          n_elts = mr->face_bi.gnum_range[1] - mr->face_bi.gnum_range[0];
        BFT_MALLOC(mr->face_gc_id, n_elts, cs_int_t);
        cs_io_read_block(&header,
                         mr->face_bi.gnum_range[0],
                         mr->face_bi.gnum_range[1],
                         mr->face_gc_id, pp_in);
      }

    }
    else if (strncmp(header.sec_name, "face_vertices_index",
                     CS_IO_NAME_LEN) == 0) {

      n_elts = mr->n_g_faces + 1;
      if (data_read != true || header.n_vals != n_elts)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else {
        fvm_gnum_t *_g_face_vertices_idx;
        cs_io_set_fvm_gnum(&header, pp_in);
        if (mr->face_bi.gnum_range[0] > 0)
          n_elts = mr->face_bi.gnum_range[1] - mr->face_bi.gnum_range[0] + 1;
        BFT_MALLOC(mr->face_vertices_idx, n_elts, fvm_lnum_t);
        BFT_MALLOC(_g_face_vertices_idx, n_elts, fvm_gnum_t);
        cs_io_read_index_block(&header,
                               mr->face_bi.gnum_range[0],
                               mr->face_bi.gnum_range[1],
                               _g_face_vertices_idx, pp_in);
        if (n_elts > 0) {
          fvm_gnum_t elt_id;
          face_vertices_idx_shift = _g_face_vertices_idx[0];
          for (elt_id = 0; elt_id < n_elts; elt_id++)
            mr->face_vertices_idx[elt_id]
              = _g_face_vertices_idx[elt_id] - face_vertices_idx_shift;
        }
        BFT_FREE(_g_face_vertices_idx);
      }

    }
    else if (strncmp(header.sec_name, "face_vertices",
                     CS_IO_NAME_LEN) == 0) {

      if (   data_read != true
          || header.n_vals != mr->n_g_face_connect_size)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else {
        size_t n_faces = mr->face_bi.gnum_range[1] - mr->face_bi.gnum_range[0];
        cs_io_set_fvm_gnum(&header, pp_in);
        n_elts =   mr->face_vertices_idx[n_faces]
                 - mr->face_vertices_idx[0];
        BFT_MALLOC(mr->face_vertices, n_elts, fvm_gnum_t);
        cs_io_read_block
          (&header,
           mr->face_vertices_idx[0] + face_vertices_idx_shift,
           mr->face_vertices_idx[n_faces] + face_vertices_idx_shift,
           mr->face_vertices, pp_in);
      }

    }
    else if (strncmp(header.sec_name, "vertex_coords",
                     CS_IO_NAME_LEN) == 0) {

      n_elts = mesh->n_g_vertices * 3;
      if (data_read != true || header.n_vals != n_elts)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else {
        if (mr->vertex_bi.gnum_range[0] > 0)
        cs_io_assert_cs_real(&header, pp_in);
          n_elts = (  mr->vertex_bi.gnum_range[1]
                    - mr->vertex_bi.gnum_range[0])*3;
        BFT_MALLOC(mr->vertex_coords, n_elts, cs_real_t);
        cs_io_read_block(&header,
                         mr->vertex_bi.gnum_range[0]*3 - 2,
                         mr->vertex_bi.gnum_range[1]*3 - 2,
                         mr->vertex_coords, pp_in);
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
        cs_io_read_global(&header, &perio_type, pp_in);
      }

    }
    else if (strncmp(header.sec_name, "periodicity_matrix_",
                     strlen("periodicity_matrix_")) == 0) {

      n_elts = 12; /* 3x4 */
      if (data_read != true || header.n_vals != n_elts)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else {
        assert(   perio_num
               == atoi(header.sec_name + strlen("periodicity_matrix_")));
        cs_io_assert_cs_real(&header, pp_in);
        cs_io_read_global(&header, perio_matrix, pp_in);

        /* Add a periodicity to mesh->periodicities */

        _add_periodicity(mesh,
                         perio_type,
                         perio_num,
                         perio_matrix);

      }

    }
    else if (strncmp(header.sec_name, "periodicity_faces_",
                     strlen("periodicity_faces_")) == 0) {

      perio_id = atoi(header.sec_name
                      + strlen("periodicity_faces_")) - 1;
      n_elts = mr->n_g_per_face_couples[perio_id] * 2;

      if (data_read != true || header.n_vals != n_elts)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name, cs_io_get_name(pp_in));
      else {

        if ((mr->per_face_bi[perio_id]).gnum_range[0] > 0)
          mr->n_per_face_couples[perio_id]
            = (  (mr->per_face_bi[perio_id]).gnum_range[1]
               - (mr->per_face_bi[perio_id]).gnum_range[0]);
        else
          mr->n_per_face_couples[perio_id] = 0;

        cs_io_set_fvm_gnum(&header, pp_in);
        n_elts = mr->n_per_face_couples[perio_id]*2;
        BFT_MALLOC(mr->per_face_couples[perio_id], n_elts, fvm_gnum_t);
        cs_io_read_block(&header,
                         (mr->per_face_bi[perio_id]).gnum_range[0]*2 -1,
                         (mr->per_face_bi[perio_id]).gnum_range[1]*2 -1,
                         mr->per_face_couples[perio_id],
                         pp_in);

      }

    }

  } /* End of loop on messages */

  /* Finalize pre-processor input */
  /*------------------------------*/

  if (cs_glob_pp_io != NULL) {
    cs_io_finalize(&cs_glob_pp_io);
    cs_glob_pp_io = NULL;
  }

  /* Read cell rank data if available */

  if (cs_glob_base_nbr > 1)
    _read_cell_rank(mesh, mr, echo);

  /* Now send data to the correct rank */
  /*-----------------------------------*/

#if defined(_CS_HAVE_MPI)

  if (cs_glob_base_nbr > 1)
    _decompose_data_g(mesh,
                      mesh_builder,
                      mr,
                      cs_glob_base_mpi_comm);

#endif

  if (cs_glob_base_nbr == 1)
    _decompose_data_l(mesh, mesh_builder, mr);

  /* Free temporary memory */

  _mesh_reader_destroy(&_cs_glob_mesh_reader);
  mr = _cs_glob_mesh_reader;
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

