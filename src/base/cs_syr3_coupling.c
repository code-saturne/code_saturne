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
 * Syrthes coupling
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#if defined(_CS_HAVE_MPI)
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

#include <fvm_nodal.h>
#include <fvm_interface.h>
#include <fvm_nodal_extract.h>
#include <fvm_nodal_project.h>
#include <fvm_nodal_triangulate.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_mesh.h"
#include "cs_mesh_connect.h"
#include "cs_mesh_select.h"
#include "cs_parall.h"
#include "cs_post.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_syr3_coupling.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

#undef _CROSS_PRODUCT_3D
#undef _MODULE_3D

enum {X, Y, Z} ;

#define _CROSS_PRODUCT_3D(cross_prod, vect1, vect2)  \
  (cross_prod[X] = vect1[Y]*vect2[Z] - vect2[Y]*vect1[Z], \
   cross_prod[Y] = vect2[X]*vect1[Z] - vect1[X]*vect2[Z], \
   cross_prod[Z] = vect1[X]*vect2[Y] - vect2[X]*vect1[Y])

#define _MODULE_3D(vect) \
  sqrt(vect[X]*vect[X] + vect[Y]*vect[Y] + vect[Z]*vect[Z])

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/* Structure associated with Syrthes coupling */

struct _cs_syr3_coupling_t {

  cs_int_t        dim;              /* Coupled mesh dimension */
  cs_int_t        ref_axis;         /* Selected axis for edge extraction */

  /* Boundary faces parameters of coupled mesh */

  cs_mesh_select_t *select;     /* Selection criteria for boundary faces */

  cs_int_t        n_coupl_faces;    /* Number of selected faces */
  cs_int_t       *coupl_face_list;  /* Local numbering of selected faces */
  cs_real_t      *weighting;        /* Triangle area or edge lengths */
  fvm_nodal_t    *coupled_mesh;     /* Nodal mesh extracted */

  fvm_interface_set_t  *if_set;     /* Interfaces between equivalent vertices */

  /* Saved arrays for post processing (float for reduced memory use) */

  cs_int_t        post_mesh_id;     /* 0 if post-processing is not active,
                                       mesh_id if post-processing is active */
  float          *wall_temp;        /* Wall temperature (received) */
  float          *flux;             /* Flux (calculated) */
  float          *tfluid_tmp;       /* Fluid temperature (points to flux in
                                       transient stage where wall_temp and
                                       fluid_temp are known, NULL once
                                       flux is calculated) */

  /* Communication structure */

  cs_syr3_comm_t      *send_comm;        /* Sending communicator */
  cs_syr3_comm_t      *recv_comm;        /* Receiving communicator */
  cs_syr3_comm_type_t  comm_type;        /* Communicator type */
  int                  comm_echo;        /* Optional echo to standard output */

#if defined (_CS_HAVE_MPI)
  cs_int_t        syr_proc_rank;    /* Syrthes rank */
#endif

};

/*============================================================================
 *  Global variables
 *============================================================================*/

static cs_int_t              cs_glob_syr_n_couplings = 0;
static cs_syr3_coupling_t  **cs_glob_syr_coupling_array = NULL;

/* Start and end (negative) numbers associated with
   dedicated post processing meshes */

static int  cs_glob_syr_post_maillage_ext[2] = {0, 1};


/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Conversion from fvm_gnum_t to cs_int_t .
 * Arrays can use the same memory.
 *
 * parameters:
 *   fvm_data            -->  input data
 *   cs_data             <--  output data
 *   n_elts              -->  numbers of elements to compute
 *----------------------------------------------------------------------------*/

static void
_convert_fvm_gnum(fvm_gnum_t   fvm_data[],
                  cs_int_t     cs_data[],
                  fvm_gnum_t   n_elts)
{
  size_t i;

  if (sizeof(cs_int_t) > sizeof(fvm_gnum_t)) {
    for (i = 0; i < n_elts; i++)
      cs_data[n_elts - 1 - i] = fvm_data[n_elts - 1 -i];
  }
  else {
    for (i = 0; i < n_elts; i++)
      cs_data[i] = fvm_data[i];
  }

}

/*----------------------------------------------------------------------------
 * Define nodal mesh for Syrthes coupling from selection criteria on border
 * faces.
 *
 * parameters:
 *   coupled_mesh_name   -->  name of the coupled mesh
 *   syr_coupling        -->  coupling structure with Syrthes
 *----------------------------------------------------------------------------*/

static fvm_nodal_t *
_define_coupled_mesh(char                *coupled_mesh_name,
                     cs_syr3_coupling_t  *syr_coupling)
{
  cs_int_t n_colors, n_groups;

  cs_mesh_select_t *selection = syr_coupling->select;
  fvm_nodal_t *coupled_mesh = NULL;

  const cs_int_t comm_echo = syr_coupling->comm_echo;

  /* Test if specific criteria exist */

  n_colors = cs_mesh_select_get_n_colors(selection);
  n_groups = cs_mesh_select_get_n_groups(selection);

  /* Creation of a new nodal mesh from selected border faces */

  if (n_groups != 0 || n_colors != 0) {

    cs_mesh_select_extract_b_faces(cs_glob_mesh,
                                   selection,
                                   &(syr_coupling->n_coupl_faces),
                                   &(syr_coupling->coupl_face_list));

    if (comm_echo >= 0) {
      bft_printf(_("\nExtracting the mesh                 ..."));
      bft_printf_flush();
    }

    coupled_mesh = cs_maillage_extrait_fac_nodal(cs_glob_mesh,
                                                 coupled_mesh_name,
                                                 0,
                                                 syr_coupling->n_coupl_faces,
                                                 NULL,
                                                 syr_coupling->coupl_face_list);

  }
  else {

    if (comm_echo >= 0) {
      bft_printf(_("\nExtracting the mesh                 ..."));
      bft_printf_flush();
    }

    coupled_mesh = cs_maillage_extrait_fac_nodal(cs_glob_mesh,
                                                 coupled_mesh_name,
                                                 0,
                                                 cs_glob_mesh->n_b_faces,
                                                 NULL,
                                                 NULL);

  }

  if (comm_echo >= 0) {
    bft_printf(" [ok]\n");
    bft_printf_flush();
  }

  return coupled_mesh;
}

/*----------------------------------------------------------------------------
 * Renumbering of selected border faces in order to be consistent with
 * parent_num from triangulation.
 *
 * parameters:
 *   syr_coupling        -->  coupling structure with Syrthes
 *----------------------------------------------------------------------------*/

static void
_renum_faces_list(cs_syr3_coupling_t *syr_coupling)
{
  cs_int_t  i, fac_id;
  cs_int_t  elt_num, elt_num_prev, n_elts;

  cs_int_t  *parent_num = NULL;
  cs_int_t  *coupl_face_list = syr_coupling->coupl_face_list;
  fvm_nodal_t  *coupled_mesh = syr_coupling->coupled_mesh;

  const int  elt_dim = syr_coupling->dim - 1;
  const cs_int_t  comm_echo = syr_coupling->comm_echo;

  if (comm_echo >= 0) {
    bft_printf(_("\n *** Renumbering of the boundary faces list    ..."));
    bft_printf_flush();
  }

  /* Retrieve parent_num from coupled mesh */

  n_elts = fvm_nodal_get_n_entities(coupled_mesh, elt_dim);

  BFT_MALLOC(parent_num, n_elts, cs_int_t);

  fvm_nodal_get_parent_num(coupled_mesh, elt_dim, parent_num);

  assert(sizeof(fvm_lnum_t) == sizeof(cs_int_t));

  /* Rebuild coupled faces list in same order as fvm_nodal_structure */

  elt_num_prev = -1;
  fac_id = 0;

  for (i = 0; i < n_elts; i++) {

    elt_num = parent_num[i];

    if (elt_num != elt_num_prev) {
      coupl_face_list[fac_id++] = elt_num;
      elt_num_prev = elt_num;
    }

  }

  assert(fac_id == syr_coupling->n_coupl_faces);

  BFT_FREE(parent_num);

  if (comm_echo >= 0) {
    bft_printf(" [ok]\n");
    bft_printf_flush();
  }

}

/*----------------------------------------------------------------------------
 * Sending to Syrthes of information about coupled vertices:
 *   - number of coupled vertices and their global numbering,
 *   - vertex coordinates.
 *
 * parameters:
 *   syr_coupling        -->  coupling structure with Syrthes
 *   n_vertices          -->  number of coupled vertices
 *   coords              <->  vertices's coordinates
 *----------------------------------------------------------------------------*/

static void
_send_coords(cs_syr3_coupling_t  *syr_coupling,
             cs_int_t             n_vertices,
             cs_real_t           *coords)
{
  cs_int_t  elt_size;

  fvm_gnum_t  n_g_vertices = 0;
  cs_int_t  _n_g_vertices = 0;
  char  *global_vtx_num_buffer = NULL;
  cs_int_t  *global_vtx_num_int = NULL;
  fvm_gnum_t  *global_vtx_num = NULL;
  fvm_nodal_t  *coupled_mesh = syr_coupling->coupled_mesh;

  const cs_int_t  dim = syr_coupling->dim;
  const cs_int_t  n_coupl_faces = syr_coupling->n_coupl_faces;

  /* Send number of vertices */

  cs_syr3_comm_envoie_message(0,
                              "coupl:b:npoinf",
                              1,
                              CS_TYPE_cs_int_t,
                              &n_vertices,
                              syr_coupling->send_comm);

  n_g_vertices = fvm_nodal_get_n_g_vertices(coupled_mesh);
  _n_g_vertices = n_g_vertices;

  /* Send global number of vertices */

  cs_syr3_comm_envoie_message(0,
                              "coupl:b:g:npoinf",
                              1,
                              CS_TYPE_cs_int_t,
                              &_n_g_vertices,
                              syr_coupling->send_comm);

  if (n_coupl_faces > 0) {

    elt_size = CS_MAX(sizeof(fvm_gnum_t), sizeof(cs_int_t));
    BFT_MALLOC(global_vtx_num_buffer, n_vertices * elt_size, char);

    global_vtx_num = (fvm_gnum_t *)global_vtx_num_buffer;

    fvm_nodal_get_global_vertex_num(coupled_mesh, global_vtx_num);

    /* Convert fvm_gnum_t to cs_int_t if necessary */

    global_vtx_num_int = (cs_int_t *)global_vtx_num_buffer;
    _convert_fvm_gnum(global_vtx_num,
                      global_vtx_num_int,
                      (fvm_gnum_t)n_vertices);

  }

  /* Send global vertex numbering */

  cs_syr3_comm_envoie_message(0,
                              "coupl:b:g:vtxnum",
                              n_vertices,
                              CS_TYPE_cs_int_t,
                              global_vtx_num_int,
                              syr_coupling->send_comm);

  if (global_vtx_num_buffer != NULL) {

    BFT_FREE(global_vtx_num_buffer);
    global_vtx_num_int = NULL;
    global_vtx_num = NULL;

  }

  /* Get vertices's coordinates */

  if (n_coupl_faces > 0) {

    /* Checkings */

    assert(sizeof(fvm_coord_t) == sizeof(cs_real_t));
    assert(sizeof(double) == sizeof(cs_real_t));

    fvm_nodal_get_vertex_coords(coupled_mesh,
                                FVM_NO_INTERLACE,
                                coords);

  }

  /* Send vertices's coordinates */

  cs_syr3_comm_envoie_message(0,
                              "coupl:b:xyzf",
                              dim * n_vertices,
                              CS_TYPE_cs_real_t,
                              coords,
                              syr_coupling->send_comm);

}

/*----------------------------------------------------------------------------
 * Send to Syrthes elements connectivity.
 *
 * parameters:
 *   syr_coupling        -->  coupling structure with Syrthes
 *   n_elts              -->  number of elements
 *----------------------------------------------------------------------------*/

static void
_send_connectivity(cs_syr3_coupling_t  *syr_coupling,
                   cs_int_t             n_elts)
{
  cs_int_t  i_elt, i_dim;
  cs_int_t  stride, elt_size;

  cs_int_t  n_connect = 0;
  char  *glob_elt_num = NULL;
  cs_int_t  *ni_connect = NULL;
  fvm_lnum_t  *connect = NULL;
  fvm_nodal_t *coupled_mesh = syr_coupling->coupled_mesh;

  const cs_int_t dim = syr_coupling->dim;
  const cs_int_t elt_dim = syr_coupling->dim - 1;
  const cs_int_t n_coupl_faces = syr_coupling->n_coupl_faces;

  /* Send number of elements */

  cs_syr3_comm_envoie_message(0,
                              "coupl:b:nelebf",
                              1,
                              CS_TYPE_cs_int_t,
                              &n_elts,
                              syr_coupling->send_comm);

  /* Get global element num */

  if (n_coupl_faces > 0) {

    elt_size = CS_MAX(sizeof(fvm_gnum_t), sizeof(cs_int_t));
    BFT_MALLOC(glob_elt_num, n_elts * elt_size, char);

    if (elt_dim == 2)
      fvm_nodal_get_global_element_num(coupled_mesh,
                                       FVM_FACE_TRIA,
                                       (fvm_gnum_t *)glob_elt_num);

    else if (elt_dim == 1)
      fvm_nodal_get_global_element_num(coupled_mesh,
                                       FVM_EDGE,
                                       (fvm_gnum_t *)glob_elt_num);

    else
      assert(elt_dim == 1 || elt_dim == 2);

    /* Convert fvm_gnum_t to cs_int_t if necessary */

    _convert_fvm_gnum((fvm_gnum_t *)glob_elt_num,
                      (cs_int_t *)glob_elt_num,
                      (fvm_gnum_t)n_elts);

  } /* n_coupl_faces > 0 */

  /* Send global element numbering */

  cs_syr3_comm_envoie_message(0,
                              "coupl:b:g:eltnum",
                              n_elts,
                              CS_TYPE_cs_int_t,
                              (cs_int_t *)glob_elt_num,
                              syr_coupling->send_comm);

  if (glob_elt_num != NULL)
    BFT_FREE(glob_elt_num);

  /* Connectivity */

  if (n_coupl_faces > 0) {

    if (elt_dim == 2) { /*If elements are triangles */

      stride = 3;
      n_connect = n_elts * stride;
      BFT_MALLOC(connect, n_connect, fvm_lnum_t);

      fvm_nodal_get_strided_connect(coupled_mesh,
                                    FVM_FACE_TRIA,
                                    connect);

    }
    else if (elt_dim == 1) { /* If elements are edges */

      stride = 2;
      n_connect = n_elts * stride;
      BFT_MALLOC(connect, n_connect, fvm_lnum_t);

      fvm_nodal_get_strided_connect(coupled_mesh,
                                    FVM_EDGE,
                                    connect);

    }
    else
      assert(elt_dim == 2 || elt_dim == 1);

    /* Convert an interlaced connectivty to a non interlaced one */

    BFT_MALLOC(ni_connect, n_connect, cs_int_t);

    for (i_dim = 0; i_dim < dim; i_dim++) {
      for (i_elt = 0; i_elt < n_elts; i_elt++)
        ni_connect[i_elt + n_elts * i_dim] = connect[i_elt * dim + i_dim];
    }

    BFT_FREE(connect);

  } /* n_coupl_faces > 0 */

  /* Send connectivity */

  cs_syr3_comm_envoie_message(0,
                              "coupl:b:nodebf",
                              n_connect,
                              CS_TYPE_cs_int_t,
                              ni_connect,
                              syr_coupling->send_comm);

  if (n_coupl_faces > 0)
    BFT_FREE(ni_connect);

}

/*----------------------------------------------------------------------------
 * Compute weighting (area for triangles and lengths for edges) used in
 * interpolation.
 *
 * parameters:
 *   syr_coupling        -->  Syrthes coupling structure
 *   coords              -->  vertices's coordinates
 *   n_elts              -->  number of elements
 *----------------------------------------------------------------------------*/

static void
_compute_weighting(cs_syr3_coupling_t  *syr_coupling,
                   cs_real_t           *coords,
                   cs_int_t             n_elts)
{
  cs_int_t  i_elt, i_stride, i_dim;
  cs_int_t  stride;

  cs_int_t  n_vertices = 0;
  cs_int_t  *connect = NULL;
  fvm_nodal_t *coupled_mesh = syr_coupling->coupled_mesh;

  const cs_int_t elt_dim = syr_coupling->dim - 1;
  const cs_int_t dim = cs_glob_mesh->dim;

  if (elt_dim == 2) { /* Triangle case */

    cs_int_t  vtx_idx[3];
    cs_real_t  vect1[3] = {0, 0, 0}, vect2[3] = {0, 0, 0};
    cs_real_t  cross_prod[3] = {0, 0, 0};

    /* Get local connectivity */

    stride = 3;
    BFT_MALLOC(connect, stride * n_elts, cs_int_t);

    fvm_nodal_get_strided_connect(coupled_mesh,
                                  FVM_FACE_TRIA,
                                  connect);

    n_vertices = (cs_int_t)fvm_nodal_get_n_entities(coupled_mesh, 0);

    /* Compute area for each triangle */

    for (i_elt = 0; i_elt < n_elts; i_elt++) {

      for (i_stride = 0; i_stride < stride; i_stride++)
        vtx_idx[i_stride] = connect[i_elt * stride + i_stride] - 1;

      /* WARNING: Syrthes coordinates are non interlaced */

      for (i_dim = 0; i_dim < dim; i_dim++) {

        vect1[i_dim] = (cs_real_t)coords[vtx_idx[1] + n_vertices * i_dim]
                     - (cs_real_t)coords[vtx_idx[0] + n_vertices * i_dim];

        vect2[i_dim] = (cs_real_t)coords[vtx_idx[2] + n_vertices * i_dim]
                     - (cs_real_t)coords[vtx_idx[0] + n_vertices * i_dim];

      }

      _CROSS_PRODUCT_3D(cross_prod, vect1, vect2);

      for (i_dim = 0; i_dim < dim; i_dim++)
        cross_prod[i_dim] *= 0.5;

      syr_coupling->weighting[i_elt] = _MODULE_3D(cross_prod);

      assert(syr_coupling->weighting[i_elt] != 0);

    } /* End of loop on elements */

  }
  else if (elt_dim == 1) { /* Edges case */

    cs_int_t  vtx_idx[2];
    cs_real_t vect[3] = {0, 0, 0};

    /* Get local connectivity */

    stride = 2;
    BFT_MALLOC(connect, stride * n_elts, cs_int_t);

    fvm_nodal_get_strided_connect(coupled_mesh,
                                  FVM_EDGE,
                                  connect);

    n_vertices = (cs_int_t)fvm_nodal_get_n_entities(coupled_mesh, 0);

    /* Compute length for each edge */

    for (i_elt = 0; i_elt < n_elts; i_elt++) {

      for (i_stride = 0; i_stride < stride; i_stride++)
        vtx_idx[i_stride] = connect[i_elt * stride + i_stride] - 1;

      /* WARNING: Syrthes coordinates are non interlaced */

      for (i_dim = 0; i_dim < dim; i_dim++)
        vect[i_dim] = (cs_real_t)coords[vtx_idx[1] + n_vertices * i_dim]
                    - (cs_real_t)coords[vtx_idx[0] + n_vertices * i_dim];

      syr_coupling->weighting[i_elt] = _MODULE_3D(vect);

      assert(syr_coupling->weighting[i_elt] != 0);

    } /* End of loop on elements */

  }
  else
    assert(elt_dim == 2 || elt_dim == 1);

  if (connect != NULL)
    BFT_FREE(connect);

}

/*----------------------------------------------------------------------------
 * Interpolate a nodal field to an element-centered field for real numbers.
 *
 * parameters:
 *   syr_coupling        -->  Syrthes coupling structure
 *   elt_values          <->  array of values defined on elements
 *   vtx_values          <->  array of values defined on vertices
 *   n_elts              -->  number of elements
 *   stride              -->  element's stride
 *   parent_num          -->  parent element numbering
 *   connect             -->  local connectivity
 *----------------------------------------------------------------------------*/

static void
_interpolate_vtx_to_elt(cs_syr3_coupling_t  *syr_coupling,
                        cs_real_t           *elt_values,
                        cs_real_t           *vtx_values,
                        cs_int_t             n_elts,
                        cs_int_t             stride,
                        cs_int_t            *parent_num,
                        cs_int_t            *connect)
{
  cs_int_t  i, j, vtx_id, fac_id;
  cs_int_t  elt_num, elt_num_prev;

  cs_real_t  *down = NULL;
  cs_real_t  up = 0;
  cs_real_t  stride_inverse = 1./stride;

  const cs_int_t  n_coupl_faces = syr_coupling->n_coupl_faces;
  const cs_real_t  *weighting = syr_coupling->weighting;

  BFT_MALLOC(down, n_coupl_faces, cs_real_t);

  /* Initialize arrays */

  for (i = 0; i < n_coupl_faces; i++) {
    elt_values[i] = 0;
    down[i] = 0;
  }

  /* Interpolate field */

  elt_num_prev = -1;
  fac_id = -1;

  for (i = 0; i < n_elts; i++) {

    elt_num = parent_num[i];

    if (elt_num != elt_num_prev) {
      fac_id += 1;
      elt_num_prev = elt_num;
    }

    up = 0.;

    for (j = 0; j < stride; j++) {
      vtx_id = connect[i*stride + j] - 1;
      up += vtx_values[vtx_id];
    }

    elt_values[fac_id] += up * stride_inverse * weighting[i];
    down[fac_id] += weighting[i];

  }

  assert(fac_id+1 == n_coupl_faces);

  for (i = 0; i < n_coupl_faces; i++)
    elt_values[i] /= down[i];

  BFT_FREE(down);

}

/*----------------------------------------------------------------------------
 * Interpolate an element-centered field to a nodal field for real numbers.
 *
 * The size of vtx_values array must be twice the number of vertices.
 * The first half gets values and the second half is used as a working array.
 * The two parts must be contiguous in parallel mode for MPI transfers.
 *
 * parameters:
 *   syr_coupling        -->  Syrthes coupling structure
 *   elt_values          <->  array of values defined on elements
 *   n_vtx_values        -->  number of values defined on vertices
 *   vtx_values          <->  array of values defined on vertices
 *   n_elts              -->  number of elements
 *   stride              -->  element's stride
 *   parent_num          -->  parent element numbering
 *   connect             -->  local connectivity
 *----------------------------------------------------------------------------*/

static void
_interpolate_elt_to_vtx(cs_syr3_coupling_t  *syr_coupling,
                        cs_real_t           *elt_values,
                        cs_int_t             n_vertices,
                        cs_real_t           *vtx_values,
                        cs_int_t             n_elts,
                        cs_int_t             stride,
                        cs_int_t            *parent_num,
                        cs_int_t            *connect)
{
  cs_int_t  i, j, fac_id, vtx_id;
  cs_int_t  elt_num, elt_num_prev;

  cs_real_t *weighting = syr_coupling->weighting;
  cs_real_t *down = vtx_values + n_vertices;

  const  cs_int_t   n_coupl_faces = syr_coupling->n_coupl_faces;

  /* Initialization of vtx_values and down => 2*n_vertices */

  for (i = 0; i < 2*n_vertices; i++)
    vtx_values[i] = 0;

  /* Compute contribution from each vertex */

  elt_num_prev = -1;
  fac_id = -1;

  for (i = 0; i < n_elts; i++) {

    elt_num = parent_num[i];

    if (elt_num != elt_num_prev) {
      fac_id += 1;
      elt_num_prev = elt_num;
    }

    for (j = 0; j < stride; j++) {

      vtx_id = connect[i*stride + j] - 1;
      vtx_values[vtx_id] += elt_values[fac_id] * weighting[i];
      down[vtx_id] += weighting[i];

    }

  }

  assert(fac_id+1 == n_coupl_faces);

  /* Sum on parallel domain boundaries ;
     Reminder: we need to have down = vtx_values + n_vtx_values */

  if (syr_coupling->if_set != NULL)
    cs_parall_interface_sr(syr_coupling->if_set, n_vertices, 2, vtx_values);

  for (i = 0; i < n_vertices; i++)
    vtx_values[i] /= down[i];

}

/*----------------------------------------------------------------------------
 * Post process variables associated with Syrthes couplings
 *
 * parameters:
 *   coupling_id         -->  Id of Syrthes coupling
 *   nt_cur_abs          -->  Current time step
 *   t_cur_abs           -->  Current time value
 *----------------------------------------------------------------------------*/

static void
_cs_syr3_coupling_post_function(cs_int_t   coupling_id,
                                cs_int_t   nt_cur_abs,
                                cs_real_t  t_cur_abs)
{
  cs_syr3_coupling_t * syr_coupling = cs_syr3_coupling_by_id(coupling_id);

  if (syr_coupling->post_mesh_id != 0) {

    cs_post_ecrit_var_som(syr_coupling->post_mesh_id,
                          _("Wall T"),
                          1,
                          false,
                          false,
                          CS_POST_TYPE_float,
                          nt_cur_abs,
                          t_cur_abs,
                          syr_coupling->wall_temp);

    cs_post_ecrit_var_som(syr_coupling->post_mesh_id,
                          _("Flux"),
                          1,
                          false,
                          false,
                          CS_POST_TYPE_float,
                          nt_cur_abs,
                          t_cur_abs,
                          syr_coupling->flux);

  }

}

/*----------------------------------------------------------------------------
 * Initialize post-processing of a Syrthes coupling
 *
 * parameters:
 *   coupling_id         -->  Id of Syrthes coupling
 *   writer_id           -->  Id of associated writer
 *----------------------------------------------------------------------------*/

static void
_cs_syr3_coupling_post_init(cs_int_t  coupling_id,
                           cs_int_t  writer_id)
{
  cs_int_t  n_vertices = 0;
  cs_int_t  mesh_id = cs_post_ret_num_maillage_libre();

  cs_syr3_coupling_t  *syr_coupling = cs_syr3_coupling_by_id(coupling_id);

  assert(syr_coupling != NULL);

  /* Exit silently if associated writer is not available */

  if (cs_post_existe_writer(writer_id) != true)
    return;

  /* Initialize post processing flag, and free previous arrays in
     case this function is called more than once */

  syr_coupling->post_mesh_id = mesh_id;

  if (syr_coupling->wall_temp != NULL)
    BFT_FREE(syr_coupling->wall_temp);

  if (syr_coupling->flux != NULL)
    BFT_FREE(syr_coupling->flux);

  /* Get number of coupled vertices */

  n_vertices = fvm_nodal_get_n_entities(syr_coupling->coupled_mesh, 0);

  /* Allocate arrays */

  if (n_vertices > 0) {
    BFT_MALLOC(syr_coupling->wall_temp, n_vertices, float);
    BFT_MALLOC(syr_coupling->flux, n_vertices, float);
  }
  syr_coupling->tfluid_tmp = NULL;

  /* Associate external mesh description with post processing subsystem */

  cs_post_ajoute_maillage_existant(mesh_id,
                                   syr_coupling->coupled_mesh,
                                   false);

  cs_post_associe(mesh_id, writer_id);

  /* Register post processing function */

  cs_post_ajoute_var_temporelle(_cs_syr3_coupling_post_function,
                                coupling_id);

  /* Update start and end (negative) numbers associated with
     dedicated post processing meshes */

  if (cs_glob_syr_post_maillage_ext[0] == 0)
    cs_glob_syr_post_maillage_ext[0] = mesh_id;

  cs_glob_syr_post_maillage_ext[1] = mesh_id;
}

/*----------------------------------------------------------------------------
 * Dump of Syrthes coupling structure
 *
 * parameters:
 *   syr_coupling        -->  Syrthes coupling structure
 *----------------------------------------------------------------------------*/

static void
_dump_syr_coupling(cs_syr3_coupling_t  *syr_coupling)
{
  int  i;

  const int  comm_echo = syr_coupling->comm_echo;

  assert(syr_coupling != NULL);

  bft_printf("\n xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
  bft_printf("   DUMP D'UNE STRUCTURE cs_syr3_syr_coupling_t\n");
  bft_printf(" xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
  bft_printf_flush();

  /* Print selection criteria */

  bft_printf("\nCritères de sélection des faces couplées\n");
  bft_printf("----------------------------------------\n");

  cs_mesh_select_dump(syr_coupling->select);

  bft_printf("\nDimension du maillage Syrthes: %i\n",
             syr_coupling->dim);

  bft_printf("Nombre de faces de bord couplées: %i\n\n",
             syr_coupling->n_coupl_faces);

  if (syr_coupling->n_coupl_faces > comm_echo) {
    for (i = 0; i < (comm_echo + 1)/2; i++)
      bft_printf("Face de bord numéro: %i\n",
                 syr_coupling->coupl_face_list[i]);
    for (i = syr_coupling->n_coupl_faces - comm_echo/2;
         i < syr_coupling->n_coupl_faces ; i++)
      bft_printf("Face de bord numéro: %i\n",
                 syr_coupling->coupl_face_list[i]);
  }
  else
    for (i = 0; i < syr_coupling->n_coupl_faces ; i++)
      bft_printf("Face de bord numéro: %i\n",
                 syr_coupling->coupl_face_list[i]);

  /* Print communicator names */

  if (syr_coupling->send_comm != NULL) {
    bft_printf("Communicateur en émission: %s\n",
               cs_syr3_comm_ret_nom(syr_coupling->send_comm));
  }
  if (syr_coupling->recv_comm != NULL) {
    bft_printf("Communicateur en réception: %s\n",
               cs_syr3_comm_ret_nom(syr_coupling->recv_comm));
  }

  bft_printf("\nType de communication: %i\n",
             syr_coupling->comm_type);

#if defined (_CS_HAVE_MPI)
  bft_printf("(MPI) rang du processus Syrthes: %i\n",
             syr_coupling->syr_proc_rank);
#endif

  bft_printf("\nxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
  bft_printf("  FIN DU DUMP DE LA STRUCTURE cs_syr3_coupling_t\n");
  bft_printf("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
  bft_printf_flush();
}

/*============================================================================
 *  Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Get number of Syrthes couplings.
 *
 * Fortran Interface:
 *
 * SUBROUTINE NBCSYR
 * *****************
 *
 * INTEGER          n_couplings     : <-- : number of Syrthes couplings
 *----------------------------------------------------------------------------*/

void CS_PROCF(nbcsyr, NBCSYR)
(
 cs_int_t  *const n_couplings
)
{
  *n_couplings = cs_glob_syr_n_couplings;
}

/*----------------------------------------------------------------------------
 * Create nodal coupled mesh.
 * Send vertices's coordinates and connectivity of coupled mesh.
 *
 * Fortran Interface:
 *
 * SUBROUTINE GEOSYR
 * *****************
 *
 * INTEGER          n_couplings     : <-- : number of Syrthes couplings
 *----------------------------------------------------------------------------*/

void CS_PROCF(geosyr, GEOSYR)
(
 cs_int_t  *const n_couplings
)
{
  cs_int_t i_coupl;

  assert(*n_couplings > 0);

  for (i_coupl = 0; i_coupl < *n_couplings; i_coupl++) {

    cs_syr3_coupling_t *syr_coupling = cs_glob_syr_coupling_array[i_coupl];

    cs_syr3_coupling_init_mesh(syr_coupling,
                               i_coupl + 1);

  }

}

/*----------------------------------------------------------------------------
 * Get number of boundary faces coupled with Syrthes.
 *
 * Fortran Interface:
 *
 * SUBROUTINE NBFSYR
 * *****************
 *
 * INTEGER          coupl_num       : --> : coupling number
 * INTEGER          n_coupl_faces   : <-- : number of coupled boundary faces
 *----------------------------------------------------------------------------*/

void CS_PROCF(nbfsyr, NBFSYR)
(
 const cs_int_t  *const coupl_num,
       cs_int_t  *const n_coupl_faces
)
{
  if (*coupl_num < 1 || *coupl_num > cs_glob_syr_n_couplings)
    bft_error(__FILE__, __LINE__, 0,
              _("SYRTHES coupling number %d impossible; "
                "there are %d couplings"),
              *coupl_num, cs_glob_syr_n_couplings);

  else
    *n_coupl_faces
      = (cs_glob_syr_coupling_array[*coupl_num - 1])->n_coupl_faces;
}


/*----------------------------------------------------------------------------
 * Get local numbering of coupled faces
 *
 * Fortran interface:
 *
 * SUBROUTINE LFASYR
 * *****************
 *
 * INTEGER          coupl_num       : --> : coupling number
 * INTEGER          n_coupl_faces   : --> : number of coupled boundary faces
 * INTEGER          coupl_face_list : <-- : list of coupled boundary faces
 *----------------------------------------------------------------------------*/

void CS_PROCF(lfasyr, LFASYR)
(
 const cs_int_t  *const coupl_num,
 const cs_int_t  *const n_coupl_faces,
       cs_int_t  *const coupl_face_list
)
{
  cs_int_t  i;
  cs_syr3_coupling_t  *coupling = NULL;

  if (*coupl_num < 1 || *coupl_num > cs_glob_syr_n_couplings)
    bft_error(__FILE__, __LINE__, 0,
              _("SYRTHES coupling number %d impossible; "
                "there are %d couplings"),
              *coupl_num, cs_glob_syr_n_couplings);
  else
    coupling = cs_glob_syr_coupling_array[*coupl_num - 1];

  if (*n_coupl_faces != coupling->n_coupl_faces)
    bft_error(__FILE__, __LINE__, 0,
              _("LFASYR: number of faces indicated for the SYRTHES\n"
                "coupling number %d incoherent (internal error)."),
              *coupl_num);

  for (i = 0 ; i < coupling->n_coupl_faces ; i++)
    coupl_face_list[i] = coupling->coupl_face_list[i];

}

/*----------------------------------------------------------------------------
 * Initialize post processing of Syrthes couplings.
 *
 * Fortran Interface:
 *
 * SUBROUTINE PSTISY
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF(pstisy, PSTISY)
(
 void
)
{
  cs_int_t i_coupl;

  for (i_coupl = 0; i_coupl < cs_glob_syr_n_couplings; i_coupl++)
    _cs_syr3_coupling_post_init(i_coupl, -1);
}

/*----------------------------------------------------------------------------
 * Get the local (negative) numbers associated with the first and last
 * post processing meshes dedicated to Syrthes couplings
 *
 * Fortran interface:
 *
 * SUBROUTINE PSTESY
 * *****************
 *
 * INTEGER          first_id        : <-- : id of first post processing mesh
 * INTEGER          last_id         : <-- : id of last post processing mesh
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstesy, PSTESY)
(
 cs_int_t  *const first_id,
 cs_int_t  *const last_id
)
{
  cs_syr3_coupling_post_id_extents(first_id, last_id);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Get number of Syrthes couplings.
 *
 * returns:
 *   number of Syrthes couplings
 *----------------------------------------------------------------------------*/

cs_int_t
cs_syr3_coupling_n_couplings(void)
{
  return cs_glob_syr_n_couplings;
}

/*----------------------------------------------------------------------------
 * Get pointer to Syrthes coupling.
 *
 * parameters:
 *   coupling_id  -->  Id (0 to n-1) of Syrthes coupling
 *
 * returns:
 *   pointer to Syrthes coupling structure
 *----------------------------------------------------------------------------*/

cs_syr3_coupling_t *
cs_syr3_coupling_by_id(cs_int_t coupling_id)
{
  cs_syr3_coupling_t  *retval = NULL;

  if (   coupling_id > -1
      && coupling_id < cs_glob_syr_n_couplings)
    retval = cs_glob_syr_coupling_array[coupling_id];

  return retval;
}

/*----------------------------------------------------------------------------
 * Get communicator type associated with Syrthes coupling
 *
 * parameters:
 *   syr_coupling        -->  Syrthes coupling structure
 *
 * returns:
 *   communicator type
 *----------------------------------------------------------------------------*/

cs_syr3_comm_type_t
cs_syr3_coupling_get_comm_type(const cs_syr3_coupling_t *syr_coupling)
{
  return syr_coupling->comm_type;
}

/*----------------------------------------------------------------------------
 * Get sending communicator associated with Syrthes coupling
 *
 * parameters:
 *   syr_coupling        -->  Syrthes coupling structure
 *
 * returns:
 *   pointer to send communicator
 *----------------------------------------------------------------------------*/

cs_syr3_comm_t *
cs_syr3_coupling_get_send_comm(const cs_syr3_coupling_t *const syr_coupling)
{
  assert(syr_coupling != NULL);

  return(syr_coupling->send_comm);
}

/*----------------------------------------------------------------------------
 * Get receiving communicator associated with Syrthes coupling
 *
 * parameters:
 *   syr_coupling        -->  Syrthes coupling structure
 *
 * returns:
 *   pointer to receive communicator
 *----------------------------------------------------------------------------*/

cs_syr3_comm_t *
cs_syr3_coupling_get_recv_comm(const cs_syr3_coupling_t *const syr_coupling)
{
  assert(syr_coupling != NULL);

  return(syr_coupling->recv_comm);
}

/*----------------------------------------------------------------------------
 * Get number of vertices in coupled mesh
 *
 * parameters:
 *   syr_coupling        -->  Syrthes coupling structure
 *
 * returns:
 *   number of vertices in coupled mesh
 *----------------------------------------------------------------------------*/

cs_int_t
cs_syr3_coupling_get_n_vertices(const cs_syr3_coupling_t *const syr_coupling)
{
  cs_int_t n_vertices = 0;

  assert(syr_coupling != NULL);

  n_vertices = fvm_nodal_get_n_entities(syr_coupling->coupled_mesh, 0);

  return n_vertices;
}

/*----------------------------------------------------------------------------
 * Create a syr_coupling_t structure
 *
 * parameters:
 *   dim                 -->  spatial mesh dimension
 *   ref_axis            -->  reference axis
 *   inv_sel             -->  invert selected faces or not
 *   n_colors            -->  number of colors
 *   colors              -->  color list
 *   n_groups            -->  number of groups
 *   groups              -->  group list
 *   syr_proc_rank       -->  syrthes processus rank
 *   comm_type           -->  communicator type
 *----------------------------------------------------------------------------*/

void
cs_syr3_coupling_add(cs_int_t       dim,
                    cs_int_t       ref_axis,
                    cs_bool_t      inv_sel,
                    cs_int_t       n_colors,
                    cs_int_t      *colors,
                    cs_int_t       n_groups,
                    char         **groups,
#if defined (_CS_HAVE_MPI)
                    cs_int_t       syr_proc_rank,
#endif
                    cs_syr3_comm_type_t comm_type)
{
  cs_syr3_coupling_t *syr_coupling = NULL;

  /* Allocate _cs_syr3_coupling_t structure */

  BFT_REALLOC(cs_glob_syr_coupling_array,
              cs_glob_syr_n_couplings + 1, cs_syr3_coupling_t*);
  BFT_MALLOC(syr_coupling, 1, cs_syr3_coupling_t);

  syr_coupling->dim = dim;
  syr_coupling->ref_axis = ref_axis;

  syr_coupling->n_coupl_faces = 0;
  syr_coupling->coupl_face_list = NULL;

  /* Selection criteria for faces */

  syr_coupling->select = cs_mesh_select_create(n_colors,
                                               n_groups,
                                               colors,
                                               groups,
                                               inv_sel);

  syr_coupling->weighting = NULL;
  syr_coupling->coupled_mesh = NULL;
  syr_coupling->if_set = NULL;

  /* Post processing */

  syr_coupling->post_mesh_id = 0;
  syr_coupling->wall_temp = NULL;
  syr_coupling->flux = NULL;

  /* Communicators */

  syr_coupling->comm_echo = -1;
  syr_coupling->comm_type = comm_type;
  syr_coupling->send_comm = NULL;
  syr_coupling->recv_comm = NULL;

#if defined (_CS_HAVE_MPI)
  syr_coupling->syr_proc_rank = syr_proc_rank;
#endif

  cs_glob_syr_coupling_array[cs_glob_syr_n_couplings] = syr_coupling;
  cs_glob_syr_n_couplings++;

}

/*----------------------------------------------------------------------------
 * Initialize communicator for Syrthes coupling
 *
 * parameters:
 *   syr_coupling        -->  Syrthes coupling structure
 *   num_syr_coupling    -->  syrthes coupling number
 *   comm_echo           -->  optional echo to standard output
 *----------------------------------------------------------------------------*/

void
cs_syr3_coupling_init_comm(cs_syr3_coupling_t *syr_coupling,
                          cs_int_t           num_syr_coupling,
                          cs_int_t           comm_echo)
{
  cs_int_t  i_coupl;

  syr_coupling->comm_echo = comm_echo;

  /* Initialize receiving communicator */

  syr_coupling->recv_comm = cs_syr3_comm_initialise("syrthes",
                                                    "solveur",
                                                    "SYRTHES_VERS_SATURNE_1.0",
                                                    num_syr_coupling,
#if defined(_CS_HAVE_MPI)
                                                    syr_coupling->syr_proc_rank,
#endif
                                                    CS_SYR3_COMM_MODE_RECEPTION,
                                                    syr_coupling->comm_type,
                                                    comm_echo);

  /* Initialize sending communicator */

  syr_coupling->send_comm = cs_syr3_comm_initialise("solveur",
                                                    "syrthes",
                                                    "SATURNE_VERS_SYRTHES_1.0",
                                                    num_syr_coupling,
#if defined(_CS_HAVE_MPI)
                                                    syr_coupling->syr_proc_rank,
#endif
                                                    CS_SYR3_COMM_MODE_EMISSION,
                                                    syr_coupling->comm_type,
                                                    comm_echo);

  if (comm_echo >= 0) {
    for (i_coupl = 0 ; i_coupl < cs_glob_syr_n_couplings; i_coupl++)
      _dump_syr_coupling(cs_glob_syr_coupling_array[i_coupl]);
  }

}

/*----------------------------------------------------------------------------
 * Destroy cs_syr3_coupling_t structures
 *----------------------------------------------------------------------------*/

void
cs_syr3_coupling_all_destroy(void)
{
  cs_int_t i_coupl;
  cs_syr3_coupling_t *syr_coupling = NULL;

  if (cs_glob_syr_n_couplings == 0)
    return;

  for (i_coupl = 0; i_coupl < cs_glob_syr_n_couplings; i_coupl++) {

    syr_coupling = cs_glob_syr_coupling_array[i_coupl];

    /* Sending "End Of File" message */

    cs_syr3_comm_envoie_message(0,
                           CS_SYR3_COMM_FIN_FICHIER,
                           0,
                           CS_TYPE_void,
                           NULL,
                           syr_coupling->send_comm);

    /* Free _cs_syr3_coupling structure */

    BFT_FREE(syr_coupling->coupl_face_list);

    /* Free post processing arrays */

    if (syr_coupling->wall_temp != NULL)
      BFT_FREE(syr_coupling->wall_temp);

    if (syr_coupling->flux != NULL)
      BFT_FREE(syr_coupling->flux);

    /* Close send and receive communicators */

    syr_coupling->send_comm = cs_syr3_comm_termine(syr_coupling->send_comm);
    syr_coupling->recv_comm = cs_syr3_comm_termine(syr_coupling->recv_comm);

    syr_coupling->select = cs_mesh_select_destroy(syr_coupling->select);

    if (syr_coupling->weighting != NULL)
      BFT_FREE(syr_coupling->weighting);

    if (syr_coupling->coupled_mesh != NULL)
      syr_coupling->coupled_mesh = fvm_nodal_destroy(syr_coupling->coupled_mesh);

    if (syr_coupling->if_set != NULL)
      syr_coupling->if_set = fvm_interface_set_destroy(syr_coupling->if_set);

#if defined(_CS_HAVE_SOCKET)
    if (syr_coupling->comm_type == CS_SYR3_COMM_TYPE_SOCKET)
      cs_syr3_comm_termine_socket();
#endif

    BFT_FREE(syr_coupling);

  } /* End of loop on cs_glob_syr_coupling_array */

  cs_glob_syr_n_couplings = 0;
  BFT_FREE(cs_glob_syr_coupling_array);

  bft_printf(_("\nDestroying the structures associated with the "
               "SYRTHES coupling.... [ok]\n"));
  bft_printf_flush();

}

/*----------------------------------------------------------------------------
 * Define coupled mesh and send it to Syrthes
 *
 * parameters:
 *   syr_coupling        -->  Syrthes coupling structure
 *   coupl_num           -->  syrthes coupling number
 *----------------------------------------------------------------------------*/

void
cs_syr3_coupling_init_mesh(cs_syr3_coupling_t *syr_coupling,
                          const cs_int_t       coupl_num)
{
  cs_int_t  length;
  cs_int_t  dim;

  fvm_gnum_t n_g_vertices = 0;
  cs_int_t   n_vertices = 0;
  cs_int_t   n_elts = 0;
  cs_int_t   n_errors = 0;

  char         *coupled_mesh_name = NULL;
  cs_real_t    *coords = NULL;
  fvm_nodal_t  *coupled_mesh = NULL;

  const cs_int_t elt_dim = syr_coupling->dim - 1;
  const cs_int_t comm_echo = syr_coupling->comm_echo;

  if (comm_echo > 0) {
    bft_printf(_("\n ** Processing the mesh for the SYRTHES coupling "
                 "\"%d\"\n\n"),
                 coupl_num);
    bft_printf_flush();
  }

  /* Define coupled mesh name */

  length = strlen("MaillageExtraitSyrthes_") + 1 + 1;
  BFT_MALLOC(coupled_mesh_name, length, char);
  sprintf(coupled_mesh_name, "MaillageExtraitSyrthes_%d", coupl_num);

  /* Define coupled mesh */

  coupled_mesh = _define_coupled_mesh(coupled_mesh_name,
                                      syr_coupling);


  n_g_vertices = fvm_nodal_get_n_g_vertices(coupled_mesh);

  if (n_g_vertices == 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Coupling with SYRTHES impossible.\n"
                "No face to couple.\n"));

  if (elt_dim == 2) {

    /* Triangulation of coupled faces */

    if (comm_echo >= 0) {
      bft_printf(_("Triangulation of the extracted mesh (%d faces)  ..."),
                 syr_coupling->n_coupl_faces);
      bft_printf_flush();
    }

    fvm_nodal_triangulate(coupled_mesh, &n_errors);

  }
  else if (elt_dim == 1) {

    /* Projection of coupled faces to edges */

    if (comm_echo >= 0) {
      bft_printf(_("Splitting the extracted mesh in edges (%d faces)  ..."),
                 syr_coupling->n_coupl_faces);
      bft_printf_flush();
    }

    fvm_nodal_project(coupled_mesh,
                      syr_coupling->ref_axis,
                      &n_errors);

  }
  else
    assert(elt_dim == 2 || elt_dim == 1);

  if (n_errors > 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Error triangulating the extracted mesh before "
                "sending to SYRTHES.\n"));

  if (comm_echo >= 0) {
    bft_printf(" [ok]\n");
    bft_printf_flush();
  }

  syr_coupling->coupled_mesh = coupled_mesh;

  /* Get number of coupled vertices */

  n_vertices  = (cs_int_t)fvm_nodal_get_n_entities(coupled_mesh, 0);

  if (syr_coupling->n_coupl_faces > 0) {

    /* Renumbering of selected border faces in order to be consistent with
       parent_num from triangulation. */

    _renum_faces_list(syr_coupling);

  } /* If n_coupl_faces > 0 */

  if (cs_glob_base_nbr > 1) {

    fvm_gnum_t  *global_vertex_num = NULL;

    BFT_MALLOC(global_vertex_num, n_vertices, fvm_gnum_t);

    fvm_nodal_get_global_vertex_num(coupled_mesh, global_vertex_num);

    /* Define interface between vertices on parallel boundaries */

    syr_coupling->if_set = fvm_interface_set_create(n_vertices,
                                                    NULL,
                                                    global_vertex_num,
                                                    NULL,
                                                    0,
                                                    NULL,
                                                    NULL,
                                                    NULL);

    BFT_FREE(global_vertex_num);

  }

  /* Communication with Syrthes */
  /*----------------------------*/

  /* Spatial dimension */

  cs_syr3_comm_envoie_message(0,
                              "coupl:b:ndim_",
                              1,
                              CS_TYPE_cs_int_t,
                              &(syr_coupling->dim),
                              syr_coupling->send_comm);

  /* Vertices information */

  dim = (cs_int_t)fvm_nodal_get_dim(coupled_mesh);

  BFT_MALLOC(coords, n_vertices * dim, cs_real_t);

  _send_coords(syr_coupling,
               n_vertices,
               coords);

  /* Element information */

  assert(elt_dim == fvm_nodal_get_max_entity_dim(coupled_mesh));
  n_elts = fvm_nodal_get_n_entities(coupled_mesh, elt_dim);

  _send_connectivity(syr_coupling,
                     n_elts);

  if (syr_coupling->n_coupl_faces > 0) {

    /*
      Compute weighting in order to interpolate element-centered values
      to vertex values (area for triangles and lengths for edges).
    */

    BFT_MALLOC(syr_coupling->weighting, n_elts, cs_real_t);

    _compute_weighting(syr_coupling,
                       coords,
                       n_elts);

  } /* n_coupl_faces > 0 */

  if (comm_echo >= 0) {

    if (elt_dim == 2)
      bft_printf(_("\nExtracted mesh built of %d triangles"),n_elts);
    else if (elt_dim == 1)
      bft_printf(_("\nExtracted mesh built of %d edges"),n_elts);
    else
      assert(elt_dim == 2 || elt_dim == 1);

    bft_printf(_(" and %d vertices (locally)\n"),n_vertices);
    bft_printf_flush();

  }

  /* Ready to start time iterations */

  cs_syr3_comm_envoie_message(0,
                              "coupl:b:debut",
                              0,
                              CS_TYPE_void,
                              NULL,
                              syr_coupling->send_comm);

  /* Free memory */

  BFT_FREE(coupled_mesh_name);

  if (coords != NULL)
    BFT_FREE(coords);

}

/*----------------------------------------------------------------------------
 * Interpolate a vertex field to an element-centered field
 *
 * parameters:
 *   syr_coupling        -->  Syrthes coupling structure
 *   vtx_values          -->  values defined on vertices
 *   elt_values          <->  values defined on elements
 *----------------------------------------------------------------------------*/

void
cs_syr3_coupling_vtx_to_elt(cs_syr3_coupling_t        *syr_coupling,
                           cs_real_t          *const vtx_values,
                           cs_real_t                *elt_values)
{
  cs_int_t n_elts, stride;

  cs_int_t  *parent_num = NULL;
  cs_int_t  *connect = NULL;

  const cs_int_t comm_echo = syr_coupling->comm_echo;
  const cs_int_t elt_dim = syr_coupling->dim - 1;
  const fvm_nodal_t *coupled_mesh = syr_coupling->coupled_mesh;

  n_elts = fvm_nodal_get_n_entities(coupled_mesh, elt_dim);

  if (n_elts == 0) return;

  /* Get parent element numbering */

  BFT_MALLOC(parent_num, n_elts, cs_int_t);
  fvm_nodal_get_parent_num(coupled_mesh, elt_dim, parent_num);

  /* Sanity test */
  assert(sizeof(fvm_lnum_t) == sizeof(cs_int_t));

  /* Get local connectivity */

  if (elt_dim == 2) { /* If elements are triangles */

    stride = 3;
    BFT_MALLOC(connect, stride * n_elts, cs_int_t);
    fvm_nodal_get_strided_connect(coupled_mesh,
                                  FVM_FACE_TRIA,
                                  connect);

  }
  else if (elt_dim == 1) { /* If elements are edges */

    stride = 2;
    BFT_MALLOC(connect, stride * n_elts, cs_int_t);
    fvm_nodal_get_strided_connect(coupled_mesh,
                                  FVM_EDGE,
                                  connect);

  }

  if (comm_echo >= 0) {
    bft_printf(_("\tInterpolation from vertices to elements            ..."));
    bft_printf_flush();
  }

  _interpolate_vtx_to_elt(syr_coupling,
                          elt_values,
                          vtx_values,
                          n_elts,
                          stride,
                          parent_num,
                          connect);

  if (comm_echo >= 0) {
    bft_printf(" [ok]\n");
    bft_printf_flush();
  }

  /* Free memory */

  BFT_FREE(connect);
  BFT_FREE(parent_num);

}

/*----------------------------------------------------------------------------
 * Interpolate an element-centered field to a vertex field.
 *
 * The size of vtx_values array must be twice the number of vertices.
 * The first half gets values and the second half is used as a working array.
 * The two parts have to be contiguous in parallel mode for MPI transfers.
 *
 * parameters:
 *   syr_coupling        -->  Syrthes coupling structure
 *   elt_values          <->  array of values defined on elements
 *   n_vtx_values        -->  number of values defined on vertices
 *   vtx_values          <->  array of values defined on vertices
 *----------------------------------------------------------------------------*/

void
cs_syr3_coupling_elt_to_vtx(cs_syr3_coupling_t        *syr_coupling,
                            cs_real_t           *const elt_values,
                            cs_int_t                   n_vertices,
                            cs_real_t                 *vtx_values)
{
  cs_int_t n_elts, stride;

  cs_int_t  *parent_num = NULL;
  cs_int_t  *connect = NULL;

  const cs_int_t elt_dim = syr_coupling->dim - 1;
  const int comm_echo = syr_coupling->comm_echo;
  const fvm_nodal_t *coupled_mesh = syr_coupling->coupled_mesh;

  /* Get number of elements */

  n_elts = fvm_nodal_get_n_entities(coupled_mesh, elt_dim);

  if (n_elts == 0) return;

  /* Get parent element numbering */

  BFT_MALLOC(parent_num, n_elts, cs_int_t);
  fvm_nodal_get_parent_num(coupled_mesh, elt_dim, parent_num);

  /* Sanity test */
  assert(sizeof(fvm_lnum_t) == sizeof(cs_int_t));

  /* Get connectivity */

  if (elt_dim == 2) { /* If elements are triangles */

    stride = 3;
    BFT_MALLOC(connect, stride * n_elts, cs_int_t);
    fvm_nodal_get_strided_connect(coupled_mesh,
                                  FVM_FACE_TRIA,
                                  connect);

  }
  else if (elt_dim == 1) { /* If elements are edges */

    stride = 2;
    BFT_MALLOC(connect, stride * n_elts, cs_int_t);
    fvm_nodal_get_strided_connect(coupled_mesh,
                                  FVM_EDGE,
                                  connect);

  }

  if (comm_echo >= 0) {
    bft_printf(_("\tInterpolation from elements to vertices            ..."));
    bft_printf_flush();
  }

  _interpolate_elt_to_vtx(syr_coupling,
                          elt_values,
                          n_vertices,
                          vtx_values,
                          n_elts,
                          stride,
                          parent_num,
                          connect);

  if (comm_echo >= 0) {
    bft_printf(" [ok]\n");
    bft_printf_flush();
  }

  /* Free memory */

  BFT_FREE(connect);
  BFT_FREE(parent_num);

}

/*----------------------------------------------------------------------------
 * Update post-processing variables of a Syrthes coupling
 *
 * parameters:
 *   syr_coupling        -->  Syrthes coupling structure
 *   step                -->  0: var = wall temperature
 *                            1: var = fluid temperature
 *                            2: var = exchange coefficient
 *   var                 -->  Pointer to variable values
 *----------------------------------------------------------------------------*/

void
cs_syr3_coupling_post_var_update(cs_syr3_coupling_t *syr_coupling,
                                int                step,
                                const cs_real_t   *var)
{
  cs_int_t  ii;
  cs_int_t  n_vertices = 0;

  assert(syr_coupling != NULL);

  if (syr_coupling->post_mesh_id == 0)
    return;

  assert(syr_coupling->wall_temp != NULL);
  assert(syr_coupling->flux != NULL);

  /* Get number of coupled vertices */

  n_vertices = fvm_nodal_get_n_entities(syr_coupling->coupled_mesh, 0);

  /* Allocate arrays */

  switch(step) {

  case 0:
    for (ii = 0; ii < n_vertices; ii++)
      syr_coupling->wall_temp[ii] = var[ii];
    break;

  case 1:
    syr_coupling->tfluid_tmp = syr_coupling->flux;
    for (ii = 0; ii < n_vertices; ii++)
      syr_coupling->tfluid_tmp[ii] = var[ii];
    break;

  case 2:
    assert(syr_coupling->tfluid_tmp == syr_coupling->flux);
    for (ii = 0; ii < n_vertices; ii++)
      syr_coupling->flux[ii] = var[ii] * (  syr_coupling->wall_temp[ii]
                                          - syr_coupling->flux[ii]);
    syr_coupling->tfluid_tmp = NULL;
    break;

  default:
    assert(0);
  }

}

/*----------------------------------------------------------------------------
 * Get the local (negative) numbers associated with the first and last
 * post processing meshes dedicated to Syrthes couplings
 *
 * parameters:
 *   first_mesh_id       <--  Id of first post processing mesh
 *   last_mesh_id        <--  Id of last post processing mesh
 *----------------------------------------------------------------------------*/

void
cs_syr3_coupling_post_id_extents(cs_int_t  *const id_mesh_start,
                                cs_int_t  *const id_mesh_end)
{
  *id_mesh_start = cs_glob_syr_post_maillage_ext[0];
  *id_mesh_end   = cs_glob_syr_post_maillage_ext[1];
}

/*----------------------------------------------------------------------------*/

/* Delete local macros */

#undef _CROSS_PRODUCT_3D
#undef _MODULE_3D

#ifdef __cplusplus
}
#endif /* __cplusplus */
