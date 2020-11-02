/*============================================================================
 * STL reader and writer.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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
#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "fvm_nodal.h"
#include "fvm_nodal_append.h"

#include "cs_mesh_headers.h"
#include "cs_order.h"
#include "cs_parall.h"
#include "cs_math.h"

/*----------------------------------------------------------------------------
 * Header of the current file
 *----------------------------------------------------------------------------*/

#include "cs_stl.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Static global variables
 *============================================================================*/

static cs_stl_mesh_info_t _stl_meshes = {NULL, 0};

cs_stl_mesh_info_t *cs_glob_stl_meshes = &_stl_meshes;

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Function to cast 4 uint8_t value into one single uint32_t value
 *
 * parameters:
 *   val <-- array of 4 uint8_t
 *
 * returns :
 *   the uint32_t corresponding value
  ----------------------------------------------------------------------------*/

static uint32_t
_cast32(uint8_t *val)
{
  uint32_t retval;

  /* We have here 4 8bits unsigned int
   * and we cast then in a single 32bit
   * unsigned int by progessive shifts */

  retval =   (uint32_t)val[0]
         + ( (uint32_t)val[1]<<8  )
         + ( (uint32_t)val[2]<<16 )
         + ( (uint32_t)val[3]<<24 );

  return retval;
}

/*----------------------------------------------------------------------------
 * Function to cut one uint32_t value into an array of 4 uint8_t value
 *
 * parameters:
 *   out --> cut array of 4 uint8_t
 *   val <-- input uint32_t value
  ----------------------------------------------------------------------------*/

static void
_cut32(uint8_t *out, int val)
{
  out[0] = (val       ) & 0xff;
  out[1] = (val >> 8  ) & 0xff;
  out[2] = (val >> 16 ) & 0xff;
  out[3] = (val >> 24 ) & 0xff;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a new STL mesh structure.
 *
 * \param[in] name name of the STL mesh
 *
 * \return pointer to the new STL mesh structure
 */
/*----------------------------------------------------------------------------*/

cs_stl_mesh_t *
cs_stl_mesh_add(const char  *name)
{
  /* First check if it already exists */
  cs_stl_mesh_t  *stl_mesh = cs_stl_mesh_get_by_name(name);

  if (stl_mesh != NULL) {
    // The STL mesh already exists
    bft_error(__FILE__, __LINE__, 0,
              _("Error creating stl mesh: mesh %s already exists."), name);

  }
  else {
    // If it does not exists create it
    _stl_meshes.n_meshes++;
    BFT_REALLOC(_stl_meshes.mesh_list, _stl_meshes.n_meshes, cs_stl_mesh_t *);

    BFT_MALLOC(stl_mesh, 1, cs_stl_mesh_t);

    if (name != NULL) {
      strncpy(stl_mesh->name, name, 9);
      stl_mesh->name[9] = '\0';
    }
    else
      bft_error(__FILE__, __LINE__, 0,
                _("Error creating stl mesh: no name given."));

    memset(stl_mesh->header, 0, 80);
    stl_mesh->n_faces = 0;
    stl_mesh->coords = NULL;
    stl_mesh->ext_mesh = NULL;

    _stl_meshes.mesh_list[_stl_meshes.n_meshes - 1] = stl_mesh;
  }

  return stl_mesh;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a pointer to a STL mesh based on its name if present.
 *
 * \param[in] name name of the STL mesh
 *
 * If no STL mesh of the given name is defined, NULL is returned.
 *
 * \return pointer to the STL mesh structure, or NULL
 */
/*----------------------------------------------------------------------------*/

cs_stl_mesh_t *
cs_stl_mesh_get_by_name(const char  *name)
{
  cs_stl_mesh_t *ptr = NULL;

  for (int s_id = 0; s_id < _stl_meshes.n_meshes; s_id ++) {
    cs_stl_mesh_t *stl_mesh = _stl_meshes.mesh_list[s_id];
    int test = strcmp(stl_mesh->name, name);
    if (test == 0)
      ptr = stl_mesh;
  }

  return ptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free all allocated STL mesh structures
 */
/*----------------------------------------------------------------------------*/

void
cs_stl_mesh_destroy_all(void)
{
  for (int s_id = 0; s_id < _stl_meshes.n_meshes; s_id ++) {
    cs_stl_mesh_t *ptr = _stl_meshes.mesh_list[s_id];
    BFT_FREE(ptr->coords);
    BFT_FREE(ptr->ext_mesh);
  }

  BFT_FREE(_stl_meshes.mesh_list);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read a binary STL file and store its content in a STL mesh structure
 *
 * Each STL file composed of the following header:
 *
 * uint8[80] – Header
 * uint32 – Number of triangles
 *
 * followed by 50 byte blocks for each triangle:
 *
 * - real32[3] – Normal vector
 * - real32[3] – Vertex 1 coordinates
 * - real32[3] – Vertex 2 coordinates
 * - real32[3] – Vertex 3 coordiantes
 * - uint16    – Attribute (any other information, usually void)
 *
 * \param[in] path     path to the STL file
 * \param[in] stl_mesh pointer to the associated STL mesh structure
 * \param[in] matrix   transformation matrix ( NULL if not used )
 */
/*----------------------------------------------------------------------------*/

void
cs_stl_file_read(cs_stl_mesh_t  *stl_mesh,
                 const char     *path,
                 double         matrix[3][4])
{
  uint8_t buf[128];
  FILE *fp;
  float *loc_coords = NULL;

  cs_lnum_t n_tria = 0;
  cs_lnum_t n_tria_new = 0;

  if (cs_glob_rank_id < 1) {

    /* Open STL file */
    fp = fopen(path, "rb");
    if (fp == NULL) {
      bft_error(__FILE__, __LINE__, 0,
                _("Error opening file \"%s\":\n\n"
                  "  %s"), path, strerror(errno));
    }

    /* Read and copy header */
    fread(buf, 84, 1, fp);
    memcpy(stl_mesh->header, buf, 80);

    /* Read number of triangles */
    uint32_t ntri;
    ntri = _cast32(buf+80);
    n_tria = (int)ntri;

    stl_mesh->n_faces = n_tria;

    /* Allocation*/
    BFT_MALLOC(stl_mesh->coords , 9*n_tria, float);
    BFT_MALLOC(loc_coords , 9, float);

    /* Loop on triangle faces
       ---------------------- */
    for (uint32_t i = 0; i < ntri; i++) {

      // Each triangle data are contained in 50bytes blocks
      fread(buf, 50, 1, fp);
      uint8_t *start = buf + 12;

      // Read the 3 vertices for the current triangle
      for (uint32_t vi = 0; vi < 3; vi ++) {

        // Read the 3 coordinates for each vertex
        for (cs_lnum_t dir = 0; dir < 3; dir ++) {

          uint32_t v_temp = _cast32(start + 3*4*vi + 4*dir);
          float *dir_coo = (float *)(&v_temp);
          loc_coords[3*vi + dir] = *dir_coo;
        }
      }

      // Check if the current triangle has a null area

      cs_real_t a[3], b[3], c[3], n[3];

      for (int dir = 0; dir < 3; dir ++) {
        a[dir] = (cs_real_t)loc_coords[3*0 + dir];
        b[dir] = (cs_real_t)loc_coords[3*1 + dir];
        c[dir] = (cs_real_t)loc_coords[3*2 + dir];
      }

      n[0] = (b[1]-a[1])*(c[2]-a[2]) - (b[2]-a[2])*(c[1]-a[1]);
      n[1] = (b[2]-a[2])*(c[0]-a[0]) - (b[0]-a[0])*(c[2]-a[2]);
      n[2] = (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]);

      cs_real_t nn = cs_math_3_norm(n);

      if (nn > 1.e-20) {
        for (int dir = 0; dir < 3; dir ++) {
          for (int vi = 0; vi < 3; vi ++) {
            stl_mesh->coords[9*n_tria_new + 3*vi + dir] = loc_coords[3*vi + dir];
          }
        }
        n_tria_new ++;
      }

    }

    /* Some log information */
    bft_printf(_("\n"
                 " ** Reading of STL mesh \"%s\"\n"
                 "    Number of triangles: %d\n"
                 "    Number of removed triangles: %d\n\n"),
               stl_mesh->name, n_tria_new, n_tria-n_tria_new);

    n_tria = n_tria_new;
    stl_mesh->n_faces = n_tria;

    /* Re-allocation*/
    BFT_REALLOC(stl_mesh->coords , 9*n_tria, float);

    BFT_FREE(loc_coords);
    fclose(fp);

    /* Apply tranformations
       -------------------- */

    if (matrix != NULL) {

      cs_real_t a[4], b[4], c[4];
      cs_real_t ar[4], br[4], cr[4];

      a[3] = 1.0;
      b[3] = 1.0;
      c[3] = 1.0;

      for (cs_lnum_t i = 0; i < n_tria; i++) {

        for (int dir = 0; dir < 4; dir ++) {
          ar[dir] = 0.0;
          br[dir] = 0.0;
          cr[dir] = 0.0;
        }

        for (int dir = 0; dir < 3; dir ++) {
          a[dir] = (cs_real_t)stl_mesh->coords[9*i + 3*0 + dir];
          b[dir] = (cs_real_t)stl_mesh->coords[9*i + 3*1 + dir];
          c[dir] = (cs_real_t)stl_mesh->coords[9*i + 3*2 + dir];
        }

        for (int dir = 0; dir < 3; dir ++) {
          for (int k = 0; k < 4; k++) {
            ar[dir] += matrix[dir][k]*a[k];
            br[dir] += matrix[dir][k]*b[k];
            cr[dir] += matrix[dir][k]*c[k];
          }
        }

        for (int dir = 0; dir < 3; dir ++) {
          stl_mesh->coords[9*i + 3*0 + dir] = (float)ar[dir];
          stl_mesh->coords[9*i + 3*1 + dir] = (float)br[dir];
          stl_mesh->coords[9*i + 3*2 + dir] = (float)cr[dir];
        }

      }
    }

  }
  /* Merge identical vertices
     ------------------------ */

  cs_lnum_t n_coo_ini = n_tria*3;

  cs_coord_t *_face_vtx_coord = NULL;
  BFT_MALLOC(_face_vtx_coord, n_coo_ini*3, cs_coord_t);
  for (cs_lnum_t i = 0; i < n_coo_ini*3; i++)
    _face_vtx_coord[i] = stl_mesh->coords[i];

  fvm_io_num_t *v_io_num = fvm_io_num_create_from_sfc(_face_vtx_coord,
                                                      3,
                                                      n_coo_ini,
                                                      FVM_IO_NUM_SFC_MORTON_BOX);

  BFT_FREE(_face_vtx_coord);

  cs_gnum_t *vtx_gnum = fvm_io_num_transfer_global_num(v_io_num);

  cs_lnum_t *order = cs_order_gnum(NULL, vtx_gnum, n_coo_ini);

  v_io_num = fvm_io_num_destroy(v_io_num);

  cs_coord_t  *vertex_coord = NULL;
  cs_lnum_t  *vertex_num = NULL;

  if (cs_glob_rank_id < 1) {

    BFT_MALLOC(vertex_coord, n_coo_ini*3, cs_coord_t);
    BFT_MALLOC(vertex_num, n_tria*3, cs_lnum_t);

    for (cs_lnum_t i = 0; i < n_coo_ini; i++)
      vertex_num[i] = -1;

    cs_lnum_t n_coo_new = 0;
    for (cs_lnum_t i = 0; i < n_coo_ini; i++) {
      cs_lnum_t j = order[i];
      vertex_num[j] = n_coo_new + 1;
      bool identical = false;
      if (! identical) {
        for (cs_lnum_t k = 0; k < 3; k++)
          vertex_coord[n_coo_new*3 + k] = stl_mesh->coords[j*3 + k];
        n_coo_new++;
      }
    }

    BFT_REALLOC(vertex_coord, n_coo_new*3, cs_coord_t);

  }

  fvm_nodal_t *ext_mesh = fvm_nodal_create(stl_mesh->name, 3);

  /* Generate FVM structure */

  assert(n_tria == 0 || cs_glob_rank_id < 1);

  fvm_nodal_append_by_transfer(ext_mesh,
                               n_tria,
                               FVM_FACE_TRIA,
                               NULL,
                               NULL,
                               NULL,
                               vertex_num,
                               NULL);

  vertex_coord = fvm_nodal_transfer_vertices(ext_mesh, vertex_coord);

  stl_mesh->ext_mesh = ext_mesh;

  /* Broadcast to other ranks */

  cs_parall_bcast(0, /* root_rank */
                  1,
                  CS_LNUM_TYPE,
                  &(stl_mesh->n_faces));

  cs_parall_bcast(0, /* root_rank */
                  stl_mesh->n_faces*9,
                  CS_FLOAT,
                  stl_mesh->coords);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Write a binary STL file according to a given STL mesh structure
 *
 * \param[in]  stl_mesh  pointer to the associated STL mesh structure
 * \param[in]  path      path to the STL file
 */
/*----------------------------------------------------------------------------*/

void
cs_stl_file_write(cs_stl_mesh_t  *stl_mesh,
                  const char     *path)
{
  uint8_t buf[128];
  FILE *fp;

  /* Output only on rank 0 for now */

  if (cs_glob_rank_id > 0)
    return;

  /* Open STL file */
  fp = fopen(path,"wb");
  if (fp == NULL) {
      bft_error(__FILE__, __LINE__, 0,
                _("Error opening file \"%s\":\n\n"
                  "  %s"), path, strerror(errno));
  }

  /* Write header */
  char header[] = "Exported from code_saturne";
  memset(buf, 0, 80);
  memcpy(buf, header, strlen(header));

  /* Cut number of triangles in 4 8bytes unsigned int */
  uint32_t ntri = (uint32_t)stl_mesh->n_faces;
  _cut32(buf+80, ntri);

  fwrite(buf, 84, 1, fp);

  /* Loop on each triangle */
  for (int i = 0; i < stl_mesh->n_faces; i++) {

    uint8_t *start = buf + 12;

    /* Compute and Write the face normal coordinates */
    float normals[3];
    float a[3], b[3], c[3];

    for (int dir = 0; dir < 3; dir ++) {
      a[dir] = stl_mesh->coords[9*i + 3*0 + dir];
      b[dir] = stl_mesh->coords[9*i + 3*1 + dir];
      c[dir] = stl_mesh->coords[9*i + 3*2 + dir];
    }

    // ab vect ac
    normals[0] = (b[1]-a[1])*(c[2]-a[2]) - (b[2]-a[2])*(c[1]-a[1]);
    normals[1] = (b[2]-a[2])*(c[0]-a[0]) - (b[0]-a[0])*(c[2]-a[2]);
    normals[2] = (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]);

    float norm = sqrt(  normals[0]*normals[0]
                      + normals[1]*normals[1]
                      + normals[2]*normals[2]);

    for (int dir = 0; dir < 3; dir ++)
      normals[dir] /= norm;

    uint32_t n_temp[3];
    for (int dir = 0; dir < 3; dir ++) {
      memcpy(&n_temp[dir],
             &normals[dir],
             sizeof n_temp[dir]);

      _cut32(buf + 4*dir, n_temp[dir]);
    }

    /* Write the 3 vertices for the current triangle */
    for (int vi = 0; vi < 3; vi ++) {

      uint32_t v_temp[3];
      for (int dir = 0; dir < 3; dir ++) {

        memcpy(&v_temp[dir],
               &stl_mesh->coords[9*i + 3*vi + dir],
               sizeof v_temp[dir]);
        _cut32(start + 3*4*vi + 4*dir, v_temp[dir]);
      }
    }

    fwrite(buf, 50, 1, fp);
  }

  fclose(fp);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

