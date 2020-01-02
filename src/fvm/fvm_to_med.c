/*============================================================================
 * Write a nodal representation associated with a mesh and associated
 * variables to MED files
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

/*----------------------------------------------------------------------------*/

#if defined(HAVE_MED)

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * MED library headers
 *----------------------------------------------------------------------------*/

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

#include <med.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"

#include "fvm_defs.h"
#include "fvm_convert_array.h"
#include "fvm_io_num.h"
#include "fvm_nodal.h"
#include "fvm_nodal_priv.h"
#include "fvm_writer_helper.h"
#include "fvm_writer_priv.h"

#include "cs_block_dist.h"
#include "cs_file.h"
#include "cs_parall.h"
#include "cs_part_to_block.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvm_to_med.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local Type Definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * MED mesh structure
 *----------------------------------------------------------------------------*/

typedef struct {

  char    name[MED_NAME_SIZE + 1];  /* Med_Mesh name */
  int     num;                      /* MED mesh number */

  med_int  entity_dim;  /* 3 for a volume mesh, 2 for a surface mesh */
  med_int  space_dim;   /* Number of coordinates to define a vertex */

} fvm_to_med_mesh_t;

/*----------------------------------------------------------------------------
 * MED field structure
 *----------------------------------------------------------------------------*/

typedef struct {

  char    name[MED_NAME_SIZE + 1];      /* MED field name */
  char    basename[MED_NAME_SIZE + 1];  /* MED field base name */

  int     id;                           /* MED field id */
  int     mesh_id;                      /* Associated mesh structure */
  int     n_components;                 /* Number of components */
  med_field_type  datatype;             /* Field datatype */

} fvm_to_med_field_t;

/*----------------------------------------------------------------------------
 * MED writer structure
 *----------------------------------------------------------------------------*/

typedef struct {

  char       *name;          /* Writer name */
  char       *filename;      /* MED file name */
  med_idt     fid;           /* MED file id */

  int                  n_med_meshes;  /* Number of MED meshes */
  fvm_to_med_mesh_t  **med_meshes;    /* Array of pointers to MED mesh
                                         structure */

  fvm_writer_time_dep_t time_dependency; /* Mesh time dependency */

  int                  n_fields;      /* Number of fields */
  fvm_to_med_field_t **fields;        /* Array of pointers to MED field
                                         structure */

  int         n_time_steps;    /* Number of meshes time steps */
  int        *time_steps;      /* Array of meshes time steps */
  double     *time_values;     /* Array of meshes time values */

  bool        allow_update;       /* Allow updates of existing data */
  bool        is_open;            /* True if MED file is open, else false */

  bool        discard_polygons;   /* Option to discard polygonal elements */
  bool        discard_polyhedra;  /* Option to discard polyhedral elements */
  bool        divide_polygons;    /* Option to tesselate polygonal elements */
  bool        divide_polyhedra;   /* Option to tesselate polyhedral elements */

  int         rank;            /* Rank of current process in communicator */
  int         n_ranks;         /* Number of processes in communicator */
  int         min_rank_step;   /* Minimum rank step for parallel IO */
  cs_lnum_t   min_block_size;  /* Minimum block size for parallel IO */

#if defined(HAVE_MPI)
  MPI_Comm    comm;            /* Associated MPI communicator */
  MPI_Comm    block_comm;      /* Associated MPI IO communicator */
#endif

} fvm_to_med_writer_t;

/*----------------------------------------------------------------------------
 * Context structure for fvm_writer_field_helper_output_* functions.
 *----------------------------------------------------------------------------*/

typedef struct {

  fvm_to_med_writer_t   *writer;        /* Pointer to writer structure */

  const char            *mesh_name;     /* med_mesh_name */
  const char            *field_name;    /* Associated MED field name */

  med_entity_type        entity_type;   /* MED entity type */
  med_geometry_type      section_type;  /* MED section type */

  int                    time_step;     /* time step number */
  double                 time_value;    /* associated time value */

  cs_gnum_t              n_g_elts;      /* Mesh global number of elements */

} _med_context_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

#define FVM_MED_MAX_N_NODES  8

static char _med_version_string_[2][32] = {"", ""};
static char _hdf5_version_string_[2][32] = {"", ""};

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Open a file associated with a writer.
 *
 * parameters:
 *   writer  <-> MED writer structure
 *   amode   <-- file access mode
 *----------------------------------------------------------------------------*/

static void
_med_file_open(fvm_to_med_writer_t  *w,
               med_access_mode       amode)
{
  assert(w->name != NULL);

  if (w->is_open)
    return;

  assert(w->fid < 0);

  if (w->allow_update)
    amode = MED_ACC_RDWR;

  /* Open in parallel when possible and useful */

#if defined(HAVE_MED_MPI)

  if (w->block_comm != MPI_COMM_NULL) {

    /* Get MPI IO hints from general code settings */

    MPI_Info hints;
    cs_file_get_default_access(CS_FILE_MODE_WRITE, NULL, &hints);

    /* Now open file */
    w->fid = MEDparFileOpen(w->filename, amode, w->block_comm, hints);
    if (w->fid < 0)
      bft_error(__FILE__, __LINE__, 0,
                _("MEDparfileOpen() failed to open file: %s"),
                w->filename);

  }
#endif

  /* Open in serial mode in other cases */

#if defined(HAVE_MPI)
  if (w->rank == 0 && w->block_comm == MPI_COMM_NULL) {
#else
  if (w->rank == 0) {
#endif
    w->fid = MEDfileOpen(w->filename, amode);
    if (w->fid < 0)
      bft_error(__FILE__, __LINE__, 0,
                _("MEDfileOpen() failed to open file: %s"),
                w->filename);
  }

  w->is_open = true;
}

/*----------------------------------------------------------------------------
 * Open a second serial only fid when needed.
 *
 * parameters:
 *   writer  <-> MED writer structure
 *----------------------------------------------------------------------------*/

static void
_med_file_open_serial(fvm_to_med_writer_t  *w)
{
  assert(w->name != NULL);

  if (w->is_open)
    return;

  assert(w->fid < 0);

  med_access_mode amode = MED_ACC_RDWR;

  /* Open in serial mode on rank 0 */

  if (w->rank == 0) {
    w->fid = MEDfileOpen(w->filename, amode);
    if (w->fid < 0)
      bft_error(__FILE__, __LINE__, 0,
                _("MEDfileOpen() failed to open file: %s"),
                w->filename);
  }

  w->is_open = true;
}

/*----------------------------------------------------------------------------
 * Close a file associated with a writer.
 *
 * parameters:
 *   writer  <-> MED writer structure
 *----------------------------------------------------------------------------*/

static void
_med_file_close(fvm_to_med_writer_t  *w)
{
  if (w->fid > -1) {
    if (MEDfileClose(w->fid) != 0)
      bft_error(__FILE__, __LINE__, 0,
                _("MEDfileClose() failed to close file \"%s\"\n"),
                w->filename);
  }

  w->fid = -1;
  w->is_open = false;
}

/*----------------------------------------------------------------------------
 * Convert FVM float type into MED float datatype according to its storage.
 *
 * parameters:
 *   fvm_data     <-- FVM data array to convert.
 *   fvm_datatype <-- FVM datatype.
 *   med_data     --> MED data array converted.
 *   n_elems      <-- Number of elements to convert.
 *
 *----------------------------------------------------------------------------*/

static void
_convert_float_fvm_to_med(const void      *fvm_data,
                          cs_datatype_t    cs_datatype,
                          med_float       *med_data,
                          int              n_elems)
{
  int  i_elem;

  const float   *data_f = (const float *)fvm_data;
  const double  *data_d = (const double *)fvm_data;

  assert(fvm_data != NULL);

  /* Type conversion adaptation */

  if (cs_datatype == CS_DOUBLE) {

    for (i_elem = 0; i_elem < n_elems; i_elem++)
      med_data[i_elem] = (med_float)data_d[i_elem];

  }
  else if (cs_datatype == CS_FLOAT)  {

    for (i_elem = 0; i_elem < n_elems; i_elem++)
      med_data[i_elem] = (med_float)data_f[i_elem];

  }
  else
    bft_error(__FILE__, __LINE__, 0,
              "_convert_float_fvm_to_med() incorrect datatype\n");

  return;
}

/*----------------------------------------------------------------------------
 * Return the MED mesh number associated with a given MED mesh name,
 * or 0 if no association found.
 *
 * parameters:
 *   writer         <-- MED writer structure
 *   med_mesh_name  <-- MED mesh name
 *
 * returns:
 *    MED mesh number, or 0 if MED mesh name is not associated with this
 *    MED writer structure in FVM
 *----------------------------------------------------------------------------*/

static int
_get_med_mesh_num(fvm_to_med_writer_t  *writer,
                  const char           *med_mesh_name)
{
  int i;
  int retval = 0;

  fvm_to_med_mesh_t  **med_meshes = writer->med_meshes;

  assert(writer != NULL);

  for (i = 0; i < writer->n_med_meshes; i++) {
    if (strcmp(med_mesh_name, med_meshes[i]->name) == 0)
      break;
  }

  if (i == writer->n_med_meshes)
    retval = 0;
  else
    retval = med_meshes[i]->num;

  return retval;
}

/*----------------------------------------------------------------------------
 * Create a MED mesh structure in MED writer structure.
 *
 * parameters:
 *   writer          <-- MED writer structure.
 *   med_mesh_name   <-- name in MED format of the mesh .
 *   mesh            <-- FVM mesh  structure.
 *
 * returns:
 *   MED mesh number, or 0 if any problem has been encountered.
 *----------------------------------------------------------------------------*/

static int
_add_med_mesh(fvm_to_med_writer_t  *writer,
              char                 *med_mesh_name,
              const fvm_nodal_t    *mesh)
{
  int i, j, id;
  int n_gc = 0;
  char med_info[MED_COMMENT_SIZE + 1] = "Generated by Code_Saturne.";
  char family_name[MED_NAME_SIZE + 1] = "";

  char  dtunit[MED_LNAME_SIZE + 1] = "s";
  char  axisname[MED_SNAME_SIZE*3 + 1];
  char  axisunit[MED_SNAME_SIZE*3 + 1];

  med_int family_num = 0;
  med_err retval = 0;

  assert(writer != NULL);

  /* Add a new MED mesh structure */

  writer->n_med_meshes += 1;
  id = writer->n_med_meshes - 1;

  BFT_REALLOC(writer->med_meshes, writer->n_med_meshes, fvm_to_med_mesh_t *);
  BFT_MALLOC(writer->med_meshes[id], 1, fvm_to_med_mesh_t);

  strncpy(writer->med_meshes[id]->name, med_mesh_name, MED_NAME_SIZE);
  writer->med_meshes[id]->name[MED_NAME_SIZE] = '\0';

  /* BUG in med: entity_dim not well treated.
     writer->med_meshes[id]->entity_dim = (med_int)entity_dim;
     So, this variable is set to space_dim */
  writer->med_meshes[id]->entity_dim = (med_int)mesh->dim;
  writer->med_meshes[id]->space_dim  = (med_int)mesh->dim;

  if (writer->fid > -1) {

    axisname[0] = 'X';
    axisname[MED_SNAME_SIZE] = 'Y';
    axisname[MED_SNAME_SIZE*2] = 'Z';
    for (j = 0; j < 3; j++)
      axisunit[MED_SNAME_SIZE*j] = 'm';
    for (i = 1; i < MED_SNAME_SIZE; i++) {
      dtunit[i] = ' ';
      for (j = 0; j < 3; j++) {
        axisname[MED_SNAME_SIZE*j + i] = ' ';
        axisunit[MED_SNAME_SIZE*j + i] = ' ';
      }
    }
    med_int space_dim = writer->med_meshes[id]->space_dim;
    dtunit[MED_LNAME_SIZE] = '\0';
    axisname[MED_SNAME_SIZE*space_dim] = '\0';
    axisunit[MED_SNAME_SIZE] = '\0';

    retval = MEDmeshCr(writer->fid,
                       med_mesh_name,
                       space_dim,
                       writer->med_meshes[id]->entity_dim,
                       MED_UNSTRUCTURED_MESH,
                       med_info,
                       dtunit,
                       MED_SORT_DTIT,
                       MED_CARTESIAN,
                       axisname,
                       axisunit);

    if (retval < 0)
      bft_error(__FILE__, __LINE__, 0,
                _("MEDmeshCr() failed to create a new med_mesh.\n"
                  "Associated med_mesh name: \"%s\"\n"
                  "Associated writer name: \"%s\"\n"),
                med_mesh_name, writer->name);

    /* Add families */

    if (mesh->gc_set != NULL)
      n_gc = fvm_group_class_set_size(mesh->gc_set);

    for (family_num = 0, retval = 0;
         - family_num < (n_gc + 1) && retval == 0;
         family_num--) {

      med_int  n_groups = 0;
      char  *med_groups = NULL;

      /* Build family name */

      snprintf(family_name, MED_NAME_SIZE, "FAMILY_%d", (int)family_num);
      family_name[MED_NAME_SIZE] = '\0';
      for (i = strlen(family_name) + 1; i < MED_NAME_SIZE; i++)
        family_name[i] = '\0';

      if (family_num < 0) {

        const fvm_group_class_t  *gc = fvm_group_class_set_get(mesh->gc_set,
                                                               - family_num - 1);
        const char **group_names = fvm_group_class_get_group_names(gc);

        /* Create MED groups */

        n_groups = fvm_group_class_get_n_groups(gc);

        if (n_groups > 0) {
          BFT_MALLOC(med_groups, MED_LNAME_SIZE  * n_groups + 1, char);
          for (i = 0; i < n_groups; i++) {
            char *group_p = med_groups + (MED_LNAME_SIZE*i);
            memset(group_p, 0, MED_LNAME_SIZE);
            strncpy(group_p, group_names[i], MED_LNAME_SIZE - 1);
          }
        }

      }

      retval = MEDfamilyCr(writer->fid,
                           med_mesh_name,
                           family_name,
                           family_num,
                           n_groups,
                           med_groups);

      if (retval < 0)
        bft_error(__FILE__, __LINE__, 0,
                  _("MEDfamilyCr() failed to create family %d"
                    "for a new med_mesh.\n"
                    "Associated med_mesh name: \"%s\"\n"
                    "Associated writer name: \"%s\"\n"),
                  (int)family_num, med_mesh_name, writer->name);

      if (med_groups != NULL)
        BFT_FREE(med_groups);

    }
  }

  writer->med_meshes[id]->num = writer->n_med_meshes;

  return writer->n_med_meshes;
}

/*----------------------------------------------------------------------------
 * Count number of extra vertices when tesselations are present
 *
 * parameters:
 *   this_writer         <-- pointer to associated writer
 *   mesh                <-- pointer to nodal mesh structure
 *   n_extra_vertices_g  --> global number of extra vertices (optional)
 *   n_extra_vertices    --> local number of extra vertices (optional)
 *----------------------------------------------------------------------------*/

static void
_count_extra_vertices(const fvm_to_med_writer_t  *this_writer,
                      const fvm_nodal_t          *mesh,
                      cs_gnum_t                  *n_extra_vertices_g,
                      cs_lnum_t                  *n_extra_vertices)
{
  int  i;

  const int  export_dim = fvm_nodal_get_max_entity_dim(mesh);

  /* Initial count and allocation */

  if (n_extra_vertices_g != NULL)
    *n_extra_vertices_g = 0;
  if (n_extra_vertices != NULL)
    *n_extra_vertices   = 0;

  for (i = 0 ; i < mesh->n_sections ; i++) {

    const fvm_nodal_section_t  *section = mesh->sections[i];

    /* Output if entity dimension equal to highest in mesh
       (i.e. no output of faces if cells present, or edges
       if cells or faces) */

    if (   section->entity_dim == export_dim
        && section->type == FVM_CELL_POLY
        && section->tesselation != NULL
        && this_writer->divide_polyhedra == true) {

      if (n_extra_vertices_g != NULL)
        *n_extra_vertices_g
          += fvm_tesselation_n_g_vertices_add(section->tesselation);

      if (n_extra_vertices != NULL)
        *n_extra_vertices
          += fvm_tesselation_n_vertices_add(section->tesselation);

    }

  }

}

/*----------------------------------------------------------------------------
 * Return extra vertex coordinates when tesselations are present
 *
 * parameters:
 *   this_writer <-- pointer to associated writer
 *   mesh        <-- pointer to nodal mesh structure that should be written
 *
 * returns:
 *   array containing all extra vertex coordinates
 *----------------------------------------------------------------------------*/

static cs_coord_t *
_extra_vertex_coords(const fvm_to_med_writer_t  *this_writer,
                     const fvm_nodal_t          *mesh)
{
  int  i;
  cs_lnum_t   n_extra_vertices_section;

  cs_lnum_t   n_extra_vertices = 0;
  size_t  coord_shift = 0;
  cs_coord_t  *coords = NULL;

  _count_extra_vertices(this_writer,
                        mesh,
                        NULL,
                        &n_extra_vertices);

  if (n_extra_vertices > 0) {

    const int  export_dim = fvm_nodal_get_max_entity_dim(mesh);

    BFT_MALLOC(coords, n_extra_vertices * 3, cs_coord_t);

    for (i = 0 ; i < mesh->n_sections ; i++) {

      const fvm_nodal_section_t  *section = mesh->sections[i];

      /* Output if entity dimension equal to highest in mesh
         (i.e. no output of faces if cells present, or edges
         if cells or faces) */

      if (   section->entity_dim == export_dim
          && section->type == FVM_CELL_POLY
          && section->tesselation != NULL
          && this_writer->divide_polyhedra == true) {

        n_extra_vertices_section
          = fvm_tesselation_n_vertices_add(section->tesselation);

        if (n_extra_vertices_section > 0) {

          fvm_tesselation_vertex_coords(section->tesselation,
                                        coords + coord_shift);

          coord_shift += n_extra_vertices_section * 3;

        }

      }
    }
  }

  return coords;
}

/*----------------------------------------------------------------------------
 * Define MED geometrical element type according to FVM element type
 *
 * parameters:
 *   fvm_elt_type  <-- pointer to fvm element type.
 *
 * return:
 *   med geometrical element type.
 *----------------------------------------------------------------------------*/

static med_geometry_type
_get_med_elt_type(const fvm_element_t fvm_elt_type)
{
  med_geometry_type  med_elt_type;

  switch (fvm_elt_type) {

  case FVM_EDGE:
    med_elt_type = MED_SEG2;
    break;

  case FVM_FACE_TRIA:
    med_elt_type = MED_TRIA3;
    break;

  case FVM_FACE_QUAD:
    med_elt_type = MED_QUAD4;
    break;

  case FVM_FACE_POLY:
    med_elt_type = MED_POLYGON;
    break;

  case FVM_CELL_TETRA:
    med_elt_type = MED_TETRA4;
    break;

  case FVM_CELL_PYRAM:
    med_elt_type = MED_PYRA5;
    break;

  case FVM_CELL_PRISM:
    med_elt_type = MED_PENTA6;
    break;

  case FVM_CELL_HEXA:
    med_elt_type = MED_HEXA8;
    break;

  case FVM_CELL_POLY:
    med_elt_type = MED_POLYHEDRON;
    break;

  default:
    med_elt_type = MED_NONE;
    bft_error(__FILE__, __LINE__, 0,
              "_get_med_elt_type(): "
              "No association with MED element type has been found\n"
              "FVM element type: \"%i\"\n",
              fvm_elt_type);

  } /* End of switch on element type */

  return med_elt_type;
}

/*----------------------------------------------------------------------------
 * Get vertex order to describe MED element type.
 *
 * parameters:
 *   med_elt_type  <-- MED element type.
 *   vertex_order  --> Pointer to vertex order array.
 *
 *----------------------------------------------------------------------------*/

static void
_get_vertex_order(const med_geometry_type med_elt_type,
                  int  *vertex_order)
{
  switch(med_elt_type) {

  case MED_SEG2:
    vertex_order[0] = 0;
    vertex_order[1] = 1;
    break;

  case MED_TRIA3:
    vertex_order[0] = 0;
    vertex_order[1] = 1;
    vertex_order[2] = 2;
    break;

  case MED_QUAD4:
    vertex_order[0] = 0;
    vertex_order[1] = 1;
    vertex_order[2] = 2;
    vertex_order[3] = 3;
    break;

  case MED_TETRA4:
    vertex_order[0] = 0;
    vertex_order[1] = 2;
    vertex_order[2] = 1;
    vertex_order[3] = 3;
    break;

  case MED_PYRA5:
    vertex_order[0] = 0;
    vertex_order[1] = 3;
    vertex_order[2] = 2;
    vertex_order[3] = 1;
    vertex_order[4] = 4;
    break;

  case MED_PENTA6:
    vertex_order[0] = 0;
    vertex_order[1] = 2;
    vertex_order[2] = 1;
    vertex_order[3] = 3;
    vertex_order[4] = 5;
    vertex_order[5] = 4;
    break;

  case MED_HEXA8:
    vertex_order[0] = 0;
    vertex_order[1] = 3;
    vertex_order[2] = 2;
    vertex_order[3] = 1;
    vertex_order[4] = 4;
    vertex_order[5] = 7;
    vertex_order[6] = 6;
    vertex_order[7] = 5;
    break;

  case MED_POLYGON:
    vertex_order[0] = -1;
    break;

  case MED_POLYHEDRON:
    vertex_order[0] = -1;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "_get_vertex_order(): No associated MED element type known\n"
              "MED element type: \"%i\"\n",
              med_elt_type);
  }

  return;
}

/*----------------------------------------------------------------------------
 * Compute the connectivity size of the current section.
 *
 * parameters:
 *   writer           <-- pointer to MED writer structure.
 *   export_sections  <-> pointer to a list of MED section structures
 *
 * returns:
 * connect_section_size <-> connectivity size of the current section.
 *----------------------------------------------------------------------------*/

static size_t
_get_connect_section_size(const fvm_to_med_writer_t   *writer,
                          const fvm_writer_section_t  *export_sections)
{
  cs_gnum_t n_g_sub_elements = 0;
  cs_gnum_t n_g_elements = 0;

  size_t connect_section_size = 0;

  const fvm_nodal_section_t *section = export_sections->section;

  if (export_sections->type == section->type) {

    /* Ordinary section */
    /*------------------*/

    if (section->stride > 0) {

      n_g_elements = fvm_nodal_section_n_g_elements(section);
      connect_section_size = n_g_elements * section->stride;

    }
    else { /* section->stride == 0 */

      cs_lnum_t   i, j, f_id;

      size_t l_connect_size = 0;

      const cs_lnum_t   *f_idx = section->face_index;
      const cs_lnum_t   *f_num = section->face_num;
      const cs_lnum_t   *v_idx = section->vertex_index;

      assert(   section->type == FVM_FACE_POLY
             || section->type == FVM_CELL_POLY);

      /* Compute local connectivity size */

      if (section->type == FVM_CELL_POLY) {

        for (i = 0; i < section->n_elements; i++) {
          for (j = f_idx[i]; j < f_idx[i+1]; j++) {
            f_id = CS_ABS(f_num[j]) - 1;
            l_connect_size += v_idx[f_id + 1] - v_idx[f_id];
          }
        }

      }
      else /* if (section->type == FVM_FACE_POLY) */
        l_connect_size = section->connectivity_size;

      if (writer->n_ranks > 1) {
#if defined(HAVE_MPI)
        MPI_Allreduce(&l_connect_size,
                      &connect_section_size,
                      1,
                      MPI_UNSIGNED_LONG,
                      MPI_SUM,
                      writer->comm);
#endif
      }
      else  /* Serial mode */
        connect_section_size = l_connect_size;

    }

  }
  else {

    /* Tesselated section */
    /*--------------------*/

    fvm_tesselation_get_global_size(export_sections->section->tesselation,
                                    export_sections->type,
                                    &n_g_sub_elements,
                                    NULL);

    connect_section_size =
      n_g_sub_elements * fvm_nodal_n_vertices_element[export_sections->type];

  }

  return connect_section_size;
}

/*----------------------------------------------------------------------------
 * Compute connectivity buffer size for MED writer file.
 *
 * parameters:
 *   writer           <-- pointer to MED writer structure.
 *   export_sections  <-> pointer to a list of MED section structures
 *
 * returns:
 *   global_connect_buffer_size: maximum buffer size useful to export
 *   connectivity.
 *----------------------------------------------------------------------------*/

static size_t
_get_connect_buffer_size(const fvm_to_med_writer_t   *writer,
                         const fvm_writer_section_t  *export_sections)

{
  med_geometry_type  med_type;

  size_t  connect_section_size = 0;
  size_t  global_connect_buffer_size = 0;

  const fvm_writer_section_t *current_section = NULL;

  med_geometry_type  previous_med_type
    = _get_med_elt_type(export_sections->type);

  current_section = export_sections;

  while (current_section != NULL) {

    med_type = _get_med_elt_type(current_section->type);

    if (med_type != previous_med_type) {

      previous_med_type = med_type;

      connect_section_size = _get_connect_section_size(writer,
                                                       current_section);

    }
    else
      connect_section_size += _get_connect_section_size(writer,
                                                        current_section);

    global_connect_buffer_size = CS_MAX(connect_section_size,
                                        global_connect_buffer_size);

    current_section = current_section->next;

  } /* End of loop on sections */

  /* When mesh has only one section */

  if (  global_connect_buffer_size == 0
     && connect_section_size > 0)
    global_connect_buffer_size = connect_section_size;

  return global_connect_buffer_size;
}

/*----------------------------------------------------------------------------
 * Adapt datatype to be MED compliant and associate its own MED datatype
 *
 * parameters:
 *   input_cs_datatype   <-- input field datatype in FVM.
 *   output_cs_datatype  --> output field datatype in FVM.
 *   med_datatype        --> associated MED field datatype.
 *   data_sizeof         --> size of datatype to be exported.
 *----------------------------------------------------------------------------*/

static void
_get_datatypes(const cs_datatype_t    input_cs_datatype,
               cs_datatype_t         *output_cs_datatype,
               med_field_type        *med_datatype,
               int                   *data_sizeof)
{
  int med_int_sizeof, med_float_sizeof;

  assert(sizeof(med_float) == 8);

  /* Verify if double define over 8 bytes */

#if (SIZEOF_DOUBLE != 8)
#error
#endif

  med_int_sizeof = sizeof(med_int);
  med_float_sizeof = sizeof(med_float);

  /* Define output_cs_datatype and med_datatype by input_cs_datatype */

  switch(input_cs_datatype) {
  case CS_DOUBLE:
    *output_cs_datatype = CS_DOUBLE;
    *med_datatype = MED_FLOAT64;
    *data_sizeof = med_float_sizeof;
    break;

  case CS_FLOAT:
    *output_cs_datatype = CS_DOUBLE;
    *med_datatype = MED_FLOAT64;
    *data_sizeof = med_float_sizeof;
    break;

  case CS_INT32:
    *output_cs_datatype = CS_INT32;
    *med_datatype = MED_INT32;
    *data_sizeof = 4;
    break;

  case CS_INT64:
    *output_cs_datatype = CS_INT32;
    *med_datatype = MED_INT32;
    *data_sizeof = 4;
    assert(0);
    break;

  case CS_UINT32:
    if (med_int_sizeof == 4) {
      *output_cs_datatype = CS_INT32;
      *med_datatype = MED_INT32;
      *data_sizeof = med_int_sizeof;
    }
    else if (med_int_sizeof == 8) {
      *output_cs_datatype = CS_INT64;
      *med_datatype = MED_INT64;
      *data_sizeof = med_int_sizeof;
    }
    else
      assert(0);
    break;

  case CS_UINT64:
    if (med_int_sizeof == 4) {
      *output_cs_datatype = CS_INT32;
      *med_datatype = MED_INT32;
      *data_sizeof = med_int_sizeof;
    }
    else if (med_int_sizeof == 8) {
      *output_cs_datatype = CS_INT64;
      *med_datatype = MED_INT64;
      *data_sizeof = med_int_sizeof;
    }
    else
      assert(0);
    break;

  default:
    assert(0);
  }

#if defined(WORKAROUND_MED_READER_INTEGER_BUG)
  /* Some versions of MEDreader crash when reading integer fields */
  *output_cs_datatype = CS_DOUBLE;
  *med_datatype = MED_FLOAT64;
#endif

  return;
}

/*----------------------------------------------------------------------------
 * Build a new field name from a base field name and a mesh name
 *
 * parameters:
 *   writer         <-- MED writer structure.
 *   med_meshname   <-- MED mesh name.
 *   base_fieldname <-- input fieldname.
 *   med_fieldname  --> MED name of the field.
 *----------------------------------------------------------------------------*/

static void
_build_new_fieldname(fvm_to_med_writer_t  *writer,
                     const char           *med_meshname,
                     const char           *base_fieldname,
                     char                 *med_fieldname)
{
  int i;
  size_t n_chars, l;

  const int n_fields = writer->n_fields - 1;

  /* Fieldname adaptation */

  strncpy(med_fieldname, base_fieldname, MED_NAME_SIZE);
  n_chars = strlen(med_fieldname);

  while (n_chars > 0) {

    if (n_chars < MED_NAME_SIZE - 4) {
      med_fieldname[n_chars] = ' ';
      med_fieldname[n_chars + 1] = '(';
    }

    strncpy(med_fieldname + n_chars + 2,
            med_meshname,
            MED_NAME_SIZE - n_chars - 3);
    med_fieldname[MED_NAME_SIZE - 1] = '\0';
    l = strlen(med_fieldname);
    med_fieldname[l] = ')';
    med_fieldname[l+1] = '\0';

    /* Loop on fields to know if field has already been created */

    for (i = 0; i < n_fields; i++) {
      fvm_to_med_field_t *field = writer->fields[i];
      if (strcmp(med_fieldname, field->name) == 0)
        break;
    }

    if (i < n_fields)  /* we have a name conflict, so start again */
      n_chars -= 1;
    else
      break;

  }

  if (n_chars < 1)
    bft_error(__FILE__, __LINE__, 0,
              _("Writer: \"%s\"\n"
                "Unable to build field name of size < %d\n"
                "for field: \"%s\" on mesh: \"%s\"."),
              writer->name, (int)MED_NAME_SIZE,
              base_fieldname, med_meshname);

  for (i = strlen(med_fieldname) + 1; i < MED_NAME_SIZE; i++)
    med_fieldname[i] = ' ';
  med_fieldname[MED_NAME_SIZE] = '\0';
}

/*----------------------------------------------------------------------------
 * Get med fieldname, update field structure and create new field if necessary
 *
 * parameters:
 *   writer          <-- MED writer structure.
 *   med_mesh_name   <-- MED mesh name.
 *   fieldname       <-- input fieldname.
 *   datatype_med    <-- associated MED field datatype.
 *   dimension       <-- dimension of field to export.
 *   med_fieldname   --> MED name of the field.
 *
 *----------------------------------------------------------------------------*/

static void
_get_med_fieldname(fvm_to_med_writer_t    *writer,
                   const char             *med_meshname,
                   const char             *fieldname,
                   med_field_type          datatype_med,
                   int                     dimension,
                   char                   *med_fieldname)
{
  int i, i_char, i_dim, n_chars, med_mesh_id;
  int n_fields, i_field, name_size;
  med_int n_components;

  char *component_name = NULL;
  char *units_name = NULL;

  char *component_unit = NULL;
  char  dt_unit[MED_LNAME_SIZE + 1] = "s";

  int basename_present = 0;

  med_err retval = 0;

  /* Fieldname adaptation */

  strncpy(med_fieldname, fieldname, MED_NAME_SIZE);
  n_chars = strlen(med_fieldname);

  for (i_char = n_chars + 1; i_char < MED_NAME_SIZE; i_char++)
    med_fieldname[i_char] = ' ';

  med_fieldname[MED_NAME_SIZE] = '\0';

  /* Get MED mesh structure */
  /*------------------------*/

  med_mesh_id = _get_med_mesh_num(writer, med_meshname) - 1;

  if (med_mesh_id < 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Writer: \"%s\"\n"
                "Mesh: \"%s\" not defined,\n"
                "but referenced by field: \"%s\"."),
              writer->name, med_meshname, med_fieldname);

  /* Loop on fields to know if field has already been created */

  n_fields = writer->n_fields;

  for (i_field = 0; i_field < n_fields; i_field++) {

    fvm_to_med_field_t *field = writer->fields[i_field];

    if (strcmp(med_fieldname, field->basename) == 0) {
      fvm_to_med_mesh_t *mesh = writer->med_meshes[field->mesh_id];
      if (strcmp(med_meshname, mesh->name) == 0)
        break;
      else
        basename_present = 1;
    }

  }

  if (i_field == n_fields) { /* Create a new field for this writer */

    BFT_REALLOC(writer->fields, writer->n_fields + 1, fvm_to_med_field_t *);

    BFT_MALLOC(writer->fields[n_fields], 1, fvm_to_med_field_t);

    memcpy(writer->fields[n_fields]->basename,
           med_fieldname,
           MED_NAME_SIZE + 1);

    if (basename_present)
      _build_new_fieldname(writer,
                           med_meshname,
                           fieldname,
                           med_fieldname); /* Updated */

    memcpy(writer->fields[n_fields]->name,
           med_fieldname,
           MED_NAME_SIZE + 1);

    writer->fields[n_fields]->id = n_fields;
    writer->fields[n_fields]->n_components = dimension;
    writer->fields[n_fields]->datatype = datatype_med;
    writer->fields[n_fields]->mesh_id = med_mesh_id;

    if (writer->fid > -1) {

      /* Component and unit names */

      name_size = MED_SNAME_SIZE * dimension;
      BFT_MALLOC(component_name, name_size + 1, char);
      BFT_MALLOC(units_name, name_size + 1, char);

      for (i = 0; i < name_size; i++) {
        component_name[i] = ' ';
        units_name[i] = ' ';
      }

      component_name[name_size] = '\0';
      units_name[name_size] = '\0';

      BFT_MALLOC(component_unit, name_size + 1, char);
      for (i = 0; i < name_size; i++)
        component_unit[i] = ' ';
      component_unit[name_size] = '\0';

      if (dimension == 1)
        sprintf(&component_name[0], "Scalar");

      else if (dimension == 3) {
        const char  *const xyz[3] = {"X", "Y", "Z"};
        for (i_dim = 0; i_dim < dimension; i_dim++)
          sprintf(&component_name[i_dim * MED_SNAME_SIZE],
                  "Component %s", xyz[i_dim]);
      }

      else if (dimension == 6) {
        const char  *const ij_sym[6] = {"11", "22", "33", "12", "13", "23"};
        for (i_dim = 0; i_dim < dimension; i_dim++)
          sprintf(&component_name[i_dim * MED_SNAME_SIZE],
                  "Component %s", ij_sym[i_dim]);
      }

      else if (dimension == 9) {
        const char  *const ij_asym[9] = {"11", "12", "13",
                                         "21", "22", "23",
                                         "31", "32", "33" };
        for (i_dim = 0; i_dim < dimension; i_dim++)
          sprintf(&component_name[i_dim * MED_SNAME_SIZE],
                  "Component %s", ij_asym[i_dim]);
      }

      for (i = 0; i < name_size; i++) {
        if (component_name[i] == '\0')
          component_name[i] = ' ';
      }

      component_name[name_size] = '\0';
      units_name[name_size] = '\0';

      /* Creation of a MED field */

      assert(writer->is_open == true);
      n_components = (med_int)dimension;

      retval = MEDfieldCr(writer->fid,
                          med_fieldname,
                          datatype_med,
                          n_components,
                          component_name,
                          component_unit,
                          dt_unit,
                          med_meshname);

      BFT_FREE(component_unit);

      if (retval < 0)
        bft_error(__FILE__, __LINE__, 0,
                  "MEDfieldCr() failed to create a field.\n"
                  "Associated writer name: \"%s\"\n"
                  "Associated mesh name: \"%s\"\n"
                  "Associated fieldname: \"%s\"\n",
                  writer->name, med_meshname, med_fieldname);

      BFT_FREE(units_name);
      BFT_FREE(component_name);

    } /* End if rank participates in IO */

    writer->n_fields++;

  } /* End of field creation */

  /* If field exists, check that dimensions and type are compatible */

  else { /*  if (i_field < n_fields) */

    fvm_to_med_field_t *field = writer->fields[i_field];

    memcpy(med_fieldname,
           writer->fields[i_field]->name,
           MED_NAME_SIZE + 1);

    if (dimension != field->n_components)
      bft_error(__FILE__, __LINE__, 0,
                _("MED field \"%s\" already defined\n"
                  "for writer \"%s\" with %d components,\n"
                  "but re-defined with %d components."),
                med_fieldname, writer->name,
                (int)field->n_components, (int)dimension);

    if (datatype_med != field->datatype)
      bft_error(__FILE__, __LINE__, 0,
                _("MED field \"%s\" already defined\n"
                  "for writer \"%s\" with datatype %d,\n"
                  "but re-defined with datatype %d."),
                med_fieldname, writer->name,
                (int)field->datatype, datatype_med);

  }

  return;
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Return datatype matching med_int
 *
 * returns:
 *   datatype matching med_float
 *----------------------------------------------------------------------------*/

static cs_datatype_t
_med_int_datatype(void)
{
  cs_datatype_t  cs_med_datatype = CS_DATATYPE_NULL;

  if (sizeof(med_int) == 4)
    cs_med_datatype = CS_INT32;
  else if (sizeof(med_int) == 8)
    cs_med_datatype = CS_INT64;
  else
    bft_error(__FILE__, __LINE__, 0 ,
              "Unexpected med_int datatype size (%d).",
              (int)(sizeof(med_float)));

  return cs_med_datatype;
}

/*----------------------------------------------------------------------------
 * Return datatype matching med_float
 *
 * returns:
 *   datatype matching med_float
 *----------------------------------------------------------------------------*/

static cs_datatype_t
_med_float_datatype(void)
{
  cs_datatype_t  cs_med_datatype = CS_DATATYPE_NULL;

  if (sizeof(med_float) == sizeof(double))
    cs_med_datatype = CS_DOUBLE;
  else if (sizeof(med_float) == sizeof(float))
    cs_med_datatype = CS_FLOAT;
  else
    bft_error(__FILE__, __LINE__, 0 ,
              "Unexpected med_float datatype size (%d).",
              (int)(sizeof(med_float)));

  return cs_med_datatype;
}

/*----------------------------------------------------------------------------
 * Write vertex coordinates to a MED file in parallel mode
 *
 * parameters:
 *   mesh          <-- pointer to nodal mesh structure that should be written.
 *   med_mesh      <-- pointer to med_mesh structure associated with the mesh.
 *   writer        <-- pointer to associated writer.
 *----------------------------------------------------------------------------*/

static void
_export_vertex_coords_g(const fvm_nodal_t    *mesh,
                        fvm_to_med_mesh_t    *med_mesh,
                        fvm_to_med_writer_t  *writer)
{
  int  i_dim;
  cs_lnum_t    i_lnod;
  cs_datatype_t  cs_datatype, cs_med_datatype;

  cs_lnum_t   n_extra_vertices = 0, n_vertices_tot = 0;
  cs_gnum_t   n_g_extra_vertices = 0, n_g_vertices_tot = 0;

  med_float  *part_coords = NULL, *block_coords = NULL;

  const int          dim = mesh->dim;
  const int          rank = writer->rank;

  cs_block_dist_info_t  bi;
  cs_part_to_block_t   *d;

  const cs_lnum_t    n_vertices = mesh->n_vertices;
  const cs_gnum_t    n_g_vertices
    = fvm_io_num_get_global_count(mesh->global_vertex_num);

  const cs_lnum_t   *parent_vertex_num = mesh->parent_vertex_num;
  const cs_coord_t  *vertex_coords = mesh->vertex_coords;

  med_err   retval = 0;

  assert(writer->is_open == true);
  assert(med_mesh != NULL);

  /* MED does not accept an empty mesh */

  if (n_g_vertices == 0)
    bft_error(__FILE__, __LINE__, 0,
              "MED does not allow to export an empty mesh,\n"
              "Mesh: \"%s\" has no vertex.\n"
              "Associated file: \"%s\".",
              mesh->name, writer->filename);

  /* Define datatypes */

  if (sizeof(cs_coord_t) == sizeof(double))
    cs_datatype = CS_DOUBLE;
  else if (sizeof(cs_coord_t) == sizeof(float))
    cs_datatype = CS_FLOAT;
  else
    bft_error(__FILE__, __LINE__, 0 ,
              "Unexpected cs_coord_t datatype size (%d).",
              (int)(sizeof(cs_coord_t)));

  cs_med_datatype = _med_float_datatype();

  /* Compute extra vertex counts if present */

  fvm_writer_count_extra_vertices(mesh,
                                  writer->divide_polyhedra,
                                  &n_g_extra_vertices,
                                  &n_extra_vertices);

  n_vertices_tot = n_vertices + n_extra_vertices;
  n_g_vertices_tot = n_g_vertices + n_g_extra_vertices;

  /* Initialize distribution info */

  fvm_writer_vertex_part_to_block_create(writer->min_rank_step,
                                         writer->min_block_size,
                                         n_g_extra_vertices,
                                         n_extra_vertices,
                                         mesh,
                                         &bi,
                                         &d,
                                         writer->comm);

  /* Build arrays */

  cs_lnum_t block_buf_size = (bi.gnum_range[1] - bi.gnum_range[0]);
  BFT_MALLOC(block_coords, block_buf_size*dim, med_float);
  BFT_MALLOC(part_coords, n_vertices_tot*dim, med_float);

  /* Export vertex coordinates to a MED file */
  /*-----------------------------------------*/

  cs_lnum_t idx = 0;

  if (parent_vertex_num != NULL) {
    for (i_lnod = 0; i_lnod < n_vertices; i_lnod++) {
      for (i_dim = 0; i_dim < dim; i_dim++)
        part_coords[idx++]
          = (med_float)vertex_coords[(parent_vertex_num[i_lnod]-1)*dim + i_dim];
    }
  }
  else
    _convert_float_fvm_to_med(vertex_coords,
                              cs_datatype,
                              part_coords,
                              n_vertices * dim);

  /* Get local extra vertex coords if present */

  {
    cs_coord_t  *extra_vertex_coords
      = fvm_writer_extra_vertex_coords(mesh, n_extra_vertices);

    idx = n_vertices*dim;

    for (i_lnod = 0; i_lnod < n_extra_vertices; i_lnod++) {
      for (i_dim = 0; i_dim < dim; i_dim++)
        part_coords[idx++] = (med_float)extra_vertex_coords[i_lnod*dim + i_dim];
    }

    BFT_FREE(extra_vertex_coords);
  }

  /* Distribute block coordinates */

  cs_part_to_block_copy_array(d,
                              cs_med_datatype,
                              dim,
                              part_coords,
                              block_coords);

  cs_part_to_block_destroy(&d);

  BFT_FREE(part_coords);

  if (writer->block_comm != MPI_COMM_NULL) { /* Parallel IO */

    med_int count = (bi.gnum_range[1] > bi.gnum_range[0]) ? 1 : 0;

    med_filter filter = MED_FILTER_INIT;
    retval = MEDfilterBlockOfEntityCr(writer->fid,
                                      n_g_vertices_tot,
                                      1,
                                      med_mesh->space_dim,
                                      MED_ALL_CONSTITUENT,
                                      MED_FULL_INTERLACE,
                                      MED_COMPACT_STMODE,
                                      MED_NO_PROFILE,
                                      bi.gnum_range[0], /* start */
                                      block_buf_size,   /* stride */
                                      count,
                                      block_buf_size,   /* blocksize */
                                      0,                /* lastblocksize */
                                      &filter);

    if (retval < 0)
      bft_error(__FILE__, __LINE__, 0,
                _("MEDfilterBlockOfEntityCr() failed for coordinates.\n"
                  "Associated writer: \"%s\"\n"
                  "Associated med_mesh: \"%s\"\n"),
                writer->name, med_mesh->name);

    retval = MEDmeshNodeCoordinateAdvancedWr(writer->fid,
                                             med_mesh->name,
                                             MED_NO_DT,
                                             MED_NO_IT,
                                             0.0,
                                             &filter,
                                             block_coords);

    if (retval < 0)
      bft_error(__FILE__, __LINE__, 0,
                _("MEDmeshNodeCoordinateAdvancedWr() failed to write coords.\n"
                  "Associated writer: \"%s\"\n"
                  "Associated med_mesh: \"%s\"\n"),
                writer->name, med_mesh->name);

    MEDfilterClose(&filter);

  }

  else if (rank == 0) {  /* Serial IO */

    if (block_coords != NULL)
      retval = MEDmeshNodeCoordinateWr(writer->fid,
                                       med_mesh->name,
                                       MED_NO_DT,
                                       MED_NO_IT,
                                       0.0,
                                       MED_FULL_INTERLACE,
                                       n_g_vertices_tot,
                                       block_coords);

    if (retval < 0)
      bft_error(__FILE__, __LINE__, 0,
                _("MEDmeshNodeCoordinateWr() failed to write coords.\n"
                  "Associated writer: \"%s\"\n"
                  "Associated med_mesh: \"%s\"\n"),
                writer->name, med_mesh->name);

  } /* End if rank == 0 */

  /* Free buffers */

  BFT_FREE(block_coords);
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Write vertex coordinates to a MED file in serial mode
 *
 * parameters:
 *   mesh          <-- pointer to nodal mesh structure that should be written.
 *   med_mesh      <-- pointer to med_mesh structure associated with the mesh.
 *   writer        <-- pointer to associated writer.
 *----------------------------------------------------------------------------*/

static void
_export_vertex_coords_l(const fvm_nodal_t     *mesh,
                        fvm_to_med_mesh_t     *med_mesh,
                        fvm_to_med_writer_t   *writer)
{
  int  i_dim;
  cs_lnum_t   i_vtx;
  cs_datatype_t  datatype;

  cs_lnum_t   idx = 0;
  cs_lnum_t   n_extra_vertices = 0;

  med_float  *med_coords = NULL;
  cs_coord_t  *extra_vertex_coords = NULL;

  const int  dim = mesh->dim;
  const cs_lnum_t   n_vertices = mesh->n_vertices;
  const cs_lnum_t    *parent_vertex_num = mesh->parent_vertex_num;
  const cs_coord_t  *vertex_coords = mesh->vertex_coords;

  med_err retval = 0;

  assert(writer->is_open == true);
  assert(med_mesh != NULL);

  /* MED does not accept an empty mesh */

  if (n_vertices == 0)
    bft_error(__FILE__, __LINE__, 0,
              "MED does not allow to export an empty mesh,\n"
              "Mesh: \"%s\" has no vertex.\n"
              "Associated file: \"%s\".",
              mesh->name, writer->filename);

  /* Define FVM datatype */

  if (sizeof(cs_coord_t) == sizeof(double))
    datatype = CS_DOUBLE;
  else if (sizeof(cs_coord_t) == sizeof(float))
    datatype = CS_FLOAT;
  else
    bft_error(__FILE__, __LINE__, 0 ,
              "Unexpected cs_coord_t datatype size (%d).",
              (int)(sizeof(cs_coord_t)));

  /* Compute extra vertex coordinates if present */

  _count_extra_vertices(writer,
                        mesh,
                        NULL,
                        &n_extra_vertices);

  extra_vertex_coords = _extra_vertex_coords(writer,
                                             mesh);

  /* Vertex coordinates export */
  /*---------------------------*/

  BFT_MALLOC(med_coords, n_vertices * dim, med_float);

  if (parent_vertex_num != NULL) {

    for (i_vtx = 0; i_vtx < n_vertices; i_vtx++) {
      for (i_dim = 0; i_dim < dim; i_dim++)
        med_coords[idx++] = (med_float)
          vertex_coords[(parent_vertex_num[i_vtx]-1) * dim + i_dim];
    }

  }
  else
    _convert_float_fvm_to_med(vertex_coords,
                              datatype,
                              med_coords,
                              n_vertices * dim);


  /* Extra vertices coordinates */
  /*----------------------------*/

  if (n_extra_vertices > 0) {

    BFT_REALLOC(med_coords,
                (n_vertices + n_extra_vertices)*dim,
                med_float);

    /* Convert cs_coord_t to med_float */

    _convert_float_fvm_to_med(extra_vertex_coords,
                              datatype,
                              med_coords + n_vertices * dim,
                              n_extra_vertices * dim);

  } /* End if n_extra_vertices > 0 */

  /* Write all coordinates */

  if (med_coords != NULL) {

    med_int _n_tot_vertices = n_vertices + n_extra_vertices;

    retval = MEDmeshNodeCoordinateWr(writer->fid,
                                     med_mesh->name,
                                     MED_NO_DT,
                                     MED_NO_IT,
                                     0.0,
                                     MED_FULL_INTERLACE,
                                     _n_tot_vertices,
                                     med_coords);

    if (retval < 0)
      bft_error(__FILE__, __LINE__, 0,
                _("MEDmeshNodeCoordinateWr() failed to write coords.\n"
                  "Associated writer: \"%s\"\n"
                  "Associated med_mesh: \"%s\"\n"),
                writer->name, med_mesh->name);

  }

  /* Free buffers */

  BFT_FREE(med_coords);

  if (extra_vertex_coords != NULL)
    BFT_FREE(extra_vertex_coords);

  return;
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Write element family values of a nodal mesh to a MED file in parallel mode.
 *
 * parameters:
 *   export_section <-- pointer to MED section helper structure
 *   writer         <-- pointer to associated writer.
 *   med_mesh       <-- pointer to MED mesh structure.
 *
 * returns:
 *  pointer to next EnSight section helper structure in list
 *----------------------------------------------------------------------------*/

static const fvm_writer_section_t *
_export_families_g(const fvm_writer_section_t  *export_section,
                   fvm_to_med_writer_t         *writer,
                   fvm_to_med_mesh_t           *med_mesh)
{
  cs_lnum_t   i, j;

  cs_block_dist_info_t  bi;
  cs_part_to_block_t  *d = NULL;

  int         n_sections = 0;
  bool        have_tesselation = false;
  cs_lnum_t   start_id = 0;
  cs_lnum_t   part_size = 0, block_size = 0;
  cs_gnum_t   block_sub_size = 0, block_end = 0;
  cs_gnum_t   n_g_elements = 0, n_g_sub_elements = 0;

  int  *part_n_sub = NULL, *block_n_sub = NULL;
  med_int  *part_values = NULL, *block_values = NULL, *_block_values = NULL;

  cs_gnum_t         *_g_elt_num = NULL;
  const cs_gnum_t   *g_elt_num
    = fvm_io_num_get_global_num(export_section->section->global_element_num);

  const fvm_writer_section_t  *current_section = NULL;

  cs_datatype_t  gc_id_type = CS_DATATYPE_NULL;
  cs_datatype_t  med_family_type = CS_DATATYPE_NULL;

  /* Initialize datatypes */

  switch(sizeof(int)) {
  case 4:
    gc_id_type = CS_INT32;
    break;
  case 8:
    gc_id_type = CS_INT64;
    break;
  default:
    assert(0);
  }

  switch(sizeof(med_int)) {
  case 4:
    med_family_type = CS_INT32;
    break;
  case 8:
    med_family_type = CS_INT64;
    break;
  default:
    assert(0);
  }

  /* Loop on sections to count output size */

  current_section = export_section;
  do {

    const fvm_nodal_section_t  *section = current_section->section;

    n_sections += 1;
    cs_gnum_t n_g_es = fvm_io_num_get_global_count(section->global_element_num);
    n_g_elements += n_g_es;
    part_size += fvm_io_num_get_local_count(section->global_element_num);
    if (current_section->type != section->type) {
      have_tesselation = true;
      cs_gnum_t n_g_ses;
      fvm_tesselation_get_global_size(section->tesselation,
                                      current_section->type,
                                      &n_g_ses,
                                      NULL);
      n_g_sub_elements += n_g_ses;
    }
    else
      n_g_sub_elements += n_g_es;

    current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous == true);

  /* Build global numbering if necessary */

  if (n_sections > 1) {

    start_id = 0;
    n_g_elements = 0;

    BFT_MALLOC(_g_elt_num, part_size, cs_gnum_t);
    g_elt_num = _g_elt_num;

    /* loop on sections which should be appended */

    current_section = export_section;
    do {

      const fvm_nodal_section_t  *section = current_section->section;
      const cs_lnum_t section_size
        = fvm_io_num_get_local_count(section->global_element_num);

      const cs_gnum_t *s_gnum
        = fvm_io_num_get_global_num(section->global_element_num);

      for (cs_lnum_t k = 0; k < section_size; k++)
        _g_elt_num[start_id + k] = s_gnum[k] + n_g_elements;

      start_id += section_size;
      n_g_elements += fvm_io_num_get_global_count(section->global_element_num);

      current_section = current_section->next;

    } while (   current_section != NULL
             && current_section->continues_previous == true);
  }

  /* Build sub-element count if necessary */

  if (have_tesselation) {

    start_id = 0;

    BFT_MALLOC(part_n_sub, part_size, int);

    current_section = export_section;
    do {

      const fvm_nodal_section_t  *section = current_section->section;
      const cs_lnum_t section_size
        = fvm_io_num_get_local_count(section->global_element_num);

      if (current_section->type != section->type) {
        const cs_lnum_t   *sub_element_idx
          = fvm_tesselation_sub_elt_index(section->tesselation,
                                          current_section->type);
        for (i = 0; i < section_size; i++)
          part_n_sub[start_id + i] = sub_element_idx[i+1] - sub_element_idx[i];
      }
      else {
        for (i = 0; i < section_size; i++)
          part_n_sub[start_id + i] = 1;
      }
      start_id += section_size;

      current_section = current_section->next;

    } while (   current_section != NULL
             && current_section->continues_previous == true);
  }

  /* Build distribution structures */

  bi = cs_block_dist_compute_sizes(writer->rank,
                                   writer->n_ranks,
                                   writer->min_rank_step,
                                   writer->min_block_size,
                                   n_g_elements);

  block_size = bi.gnum_range[1] - bi.gnum_range[0];

  d = cs_part_to_block_create_by_gnum(writer->comm, bi, part_size, g_elt_num);

  if (_g_elt_num != NULL)
    cs_part_to_block_transfer_gnum(d, _g_elt_num);

  g_elt_num = NULL;
  _g_elt_num = NULL;

  /* Distribute sub-element info in case of tesselation */

  if (have_tesselation) {

    BFT_MALLOC(block_n_sub, block_size, int);
    cs_part_to_block_copy_array(d,
                                CS_INT_TYPE,
                                1,
                                part_n_sub,
                                block_n_sub);
    BFT_FREE(part_n_sub);

    for (i = 0; i < block_size; i++)
      block_sub_size += block_n_sub[i];

  }
  else
    block_sub_size = block_size;

  /* To save space, in case of tesselation, part_values and _block_values
     point to the same memory space, as they are not needed simultaneously.
     Without tesselation, _block_values simply points to block_values */

  BFT_MALLOC(block_values, block_size, med_int);

  if (have_tesselation) {
    BFT_MALLOC(part_values,
               CS_MAX(part_size, (cs_lnum_t)block_sub_size),
               med_int);
    MPI_Scan(&block_sub_size, &block_end, 1, CS_MPI_GNUM, MPI_SUM,
             writer->comm);
    block_end += 1;
    _block_values = part_values;
  }
  else {
    BFT_MALLOC(part_values, part_size, med_int);
    block_end = bi.gnum_range[1];
    _block_values = block_values;
  }

  /* Distribute partition to block values */

  /* loop on sections which should be appended */

  start_id = 0;

  current_section = export_section;
  do {

    const fvm_nodal_section_t  *section = current_section->section;
    const void *src_data[1] = {section->gc_id};

    if (section->gc_id != NULL) {

      fvm_convert_array(1, /* src_dim */
                        0, /* src_dim_shift */
                        1, /* dest_dim */
                        0,
                        section->n_elements,
                        CS_NO_INTERLACE,
                        gc_id_type,
                        med_family_type,
                        0, /* n_parent_lists */
                        NULL,
                        NULL,
                        src_data,
                        part_values + start_id);

      for (i = 0; i < section->n_elements; i++)
        part_values[start_id + i] *= -1;

    }
    else {

      for (i = 0; i < section->n_elements; i++)
        part_values[start_id + i] = 0;

    }

    start_id += fvm_io_num_get_local_count(section->global_element_num);

    current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous == true);

  /* Distribute part values */

  cs_part_to_block_copy_array(d,
                              med_family_type,
                              1,
                              part_values,
                              block_values);

  /* Scatter values to sub-elements in case of tesselation */

  if (have_tesselation) {
    cs_lnum_t   k = 0;
    for (i = 0; i < block_size; i++) {
      for (j = 0; j < block_n_sub[i]; j++)
        _block_values[k++] = block_values[i];
    }
  }

  med_geometry_type  med_section_type
    = _get_med_elt_type(export_section->type);

  /* Write block values */

  if (writer->block_comm != MPI_COMM_NULL) { /* Parallel IO */

    med_err retval;
    med_int count = (bi.gnum_range[1] > bi.gnum_range[0]) ? 1 : 0;
    cs_gnum_t block_start = block_end - block_sub_size;

    med_filter filter = MED_FILTER_INIT;
    retval = MEDfilterBlockOfEntityCr(writer->fid,
                                      n_g_sub_elements,
                                      1,
                                      1,
                                      MED_ALL_CONSTITUENT,
                                      MED_FULL_INTERLACE,
                                      MED_COMPACT_STMODE,
                                      MED_NO_PROFILE,
                                      block_start,      /* start */
                                      block_sub_size,   /* stride */
                                      count,
                                      block_sub_size,   /* blocksize */
                                      0,                /* lastblocksize */
                                      &filter);

    if (retval < 0)
      bft_error(__FILE__, __LINE__, 0,
                _("MEDfilterBlockOfEntityCr() failed for families.\n"
                  "Associated writer: \"%s\"\n"
                  "Associated med_mesh: \"%s\"\n"),
                writer->name, med_mesh->name);

    retval = MEDmeshEntityAttributeAdvancedWr(writer->fid,
                                              med_mesh->name,
                                              MED_FAMILY_NUMBER,
                                              MED_NO_DT,
                                              MED_NO_IT,
                                              MED_CELL,
                                              med_section_type,
                                              &filter,
                                              _block_values);

    if (retval < 0)
      bft_error
        (__FILE__, __LINE__, 0,
         _("MEDmeshEntityAttributeAdvancedWr() failed to write family numbers:\n"
           "Associated writer: \"%s\"\n"
           "Associated med_mesh_name: \"%s\"\n"
           "Associated MED geometrical element: \"%i\"\n"),
         writer->name, med_mesh->name, med_section_type);

    MEDfilterClose(&filter);

  }
  else if (writer->rank == 0) {

    med_err retval = MEDmeshEntityFamilyNumberWr(writer->fid,
                                                 med_mesh->name,
                                                 MED_NO_DT,
                                                 MED_NO_IT,
                                                 MED_CELL,
                                                 med_section_type,
                                                 block_sub_size,
                                                 _block_values);

    if (retval < 0)
      bft_error
        (__FILE__, __LINE__, 0,
         _("MEDmeshEntityFamilyNumberWr() failed to write family numbers:\n"
           "Associated writer: \"%s\"\n"
           "Associated med_mesh_name: \"%s\"\n"
           "Associated MED geometrical element: \"%i\"\n"),
         writer->name, med_mesh->name, med_section_type);

  }

  BFT_FREE(block_values);
  BFT_FREE(part_values);

  cs_part_to_block_destroy(&d);

  if (block_n_sub != NULL)
    BFT_FREE(block_n_sub);

  return current_section;
}

#endif /* HAVE_MPI */

/*----------------------------------------------------------------------------
 * Count local and global elements for a given section type.
 *
 * parameters:
 *   export_sections <-- pointer to sections list to export.
 *   n_elements      --> number of elements for this section type
 *   n_g_elements    --> global number of elements for this section type
 *----------------------------------------------------------------------------*/

static void
_count_connect_g(const fvm_writer_section_t  *export_sections,
                 cs_lnum_t                   *n_elements,
                 cs_gnum_t                   *n_g_elements)
{
  const fvm_writer_section_t *current_section = export_sections;

  *n_elements = 0;
  *n_g_elements = 0;

  /* Compute cumulative sizes of sections sharing the same element type */
  /*--------------------------------------------------------------------*/

  do {   /* Loop on sections with equivalent MED element type */

    const fvm_nodal_section_t *section = current_section->section;

    if (section->type == current_section->type) { /* Normal section */
      *n_elements += section->n_elements;
      *n_g_elements += fvm_nodal_section_n_g_elements(section);
    }
    else { /* Tesselated section */
      cs_gnum_t n_g_sub_elements = 0;
      *n_elements += fvm_tesselation_n_sub_elements(section->tesselation,
                                                    current_section->type);
      fvm_tesselation_get_global_size(section->tesselation,
                                      current_section->type,
                                      &n_g_sub_elements,
                                      NULL);
      *n_g_elements += n_g_sub_elements;
    }

    current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous);
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Return global element number array if needed.
 *
 * If there is only one section of the current element type, and it is not a
 * tesselated section, the returned array will be NULL, as a pointer to
 * fvm_io_num_get_global_num(export_sections->section->global_element_num)
 * is enough.
 *
 * The caller is responsible for freeing the returned array.
 *
 * parameters:
 *   export_sections       <-- pointer to sections list to export.
 *
 * returns:
 *   pointer to section or allocated element global numbers
 *----------------------------------------------------------------------------*/

static cs_gnum_t *
_section_elt_gnum(const fvm_writer_section_t  *export_sections)
{
  bool have_tesselation = false;
  cs_lnum_t n_elements = 0;
  const fvm_writer_section_t *current_section = NULL;

  cs_gnum_t *elt_gnum = NULL;

  /* Compute cumulative sizes of sections sharing the same element type */

  current_section = export_sections;

  do {   /* Loop on sections with equivalent MED element type */

    const fvm_nodal_section_t *section = current_section->section;

    if (section->type == current_section->type)
      n_elements += section->n_elements;
    else {
      n_elements += fvm_tesselation_n_sub_elements(section->tesselation,
                                                   current_section->type);
      have_tesselation = true;
    }

    current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous);

  /* Single section with no tesselation case */

  if (!have_tesselation && export_sections->section->n_elements == n_elements)
    return elt_gnum;

  /* Case where the array must be assembled */

  BFT_MALLOC(elt_gnum, n_elements, cs_gnum_t);

  cs_lnum_t elt_id = 0;
  cs_gnum_t elt_gnum_shift = 0;

  current_section = export_sections;

  do {   /* Loop on sections with equivalent MED element type */

    const fvm_nodal_section_t *section = current_section->section;

    if (section->type == current_section->type) {
      const cs_gnum_t *s_elt_gnum
        = fvm_io_num_get_global_num(section->global_element_num);
      for (cs_lnum_t i = 0; i < section->n_elements; i++)
        elt_gnum[elt_id++] = s_elt_gnum[i] + elt_gnum_shift;
      elt_gnum_shift += fvm_io_num_get_global_count(section->global_element_num);
    }
    else {
      cs_lnum_t n_s_elements
        = fvm_tesselation_n_sub_elements(section->tesselation,
                                         current_section->type);
      const cs_lnum_t *sub_index
        = fvm_tesselation_sub_elt_index(section->tesselation,
                                        current_section->type);
      cs_lnum_t *n_sub_entities;
      BFT_MALLOC(n_sub_entities, section->n_elements, cs_lnum_t);
      for (cs_lnum_t i = 0; i < section->n_elements; i++)
        n_sub_entities[i] = sub_index[i+1] - sub_index[i];
      fvm_io_num_t *sub_io_num
        = fvm_io_num_create_from_sub(section->global_element_num,
                                     n_sub_entities);
      BFT_FREE(n_sub_entities);
      const cs_gnum_t *s_elt_gnum
        = fvm_io_num_get_global_num(sub_io_num);
      for (cs_lnum_t i = 0; i < n_s_elements; i++)
        elt_gnum[elt_id++] = s_elt_gnum[i] + elt_gnum_shift;
      elt_gnum_shift += fvm_io_num_get_global_count(sub_io_num);
      sub_io_num = fvm_io_num_destroy(sub_io_num);
    }

    current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous);

  return elt_gnum;
}

/*----------------------------------------------------------------------------
 * Write strided elements connectivity to a MED file in parallel mode
 *
 * parameters:
 *   export_sections       <-- pointer to sections list to export.
 *   writer                <-- pointer to associated writer.
 *   mesh                  <-- pointer to FVM mesh structure.
 *   med_mesh              <-- pointer to MED mesh structure.
 *
 * returns:
 *  pointer to next MED section structure in list
 *----------------------------------------------------------------------------*/

static const fvm_writer_section_t *
_export_connect_g(const fvm_writer_section_t  *export_sections,
                  fvm_to_med_writer_t         *w,
                  const fvm_nodal_t           *mesh,
                  fvm_to_med_mesh_t           *med_mesh)
{
  int vertex_order[FVM_MED_MAX_N_NODES];
  cs_lnum_t   l_id, elt_id;
  med_geometry_type  med_section_type;

  const int stride = fvm_nodal_n_vertices_element[export_sections->type];
  const fvm_writer_section_t *current_section = NULL;
  const cs_gnum_t *g_vtx_num
    = fvm_io_num_get_global_num(mesh->global_vertex_num);

  med_int *part_vertex_num = NULL, *block_vertex_num = NULL;

  med_err  retval = 0;

  current_section = export_sections;

  /* Get MED element vertex order */

  med_section_type = _get_med_elt_type(current_section->type);

  _get_vertex_order(med_section_type, vertex_order);

  /* Count elements and prepare block distribution */

  cs_lnum_t   n_part_elts = 0;
  cs_gnum_t   n_g_elts = 0;

  _count_connect_g(export_sections, &n_part_elts, &n_g_elts);

  const cs_block_dist_info_t bi = cs_block_dist_compute_sizes(w->rank,
                                                              w->n_ranks,
                                                              w->min_rank_step,
                                                              w->min_block_size,
                                                              n_g_elts);

  const cs_lnum_t n_block_elts = bi.gnum_range[1] - bi.gnum_range[0];

  const cs_gnum_t *s_elt_gnum
    = fvm_io_num_get_global_num(export_sections->section->global_element_num);
  cs_gnum_t *_s_elt_gnum = _section_elt_gnum(export_sections);

  if (_s_elt_gnum != NULL)
    s_elt_gnum = _s_elt_gnum;

  /* Distribute connectivity from sections sharing the same element type */
  /*---------------------------------------------------------------------*/

  BFT_MALLOC(block_vertex_num, stride * n_block_elts, med_int);
  BFT_MALLOC(part_vertex_num, stride * n_part_elts, med_int);

  cs_part_to_block_t *d
    = cs_part_to_block_create_by_gnum(w->comm, bi, n_part_elts, s_elt_gnum);

  if (_s_elt_gnum != NULL)
    cs_part_to_block_transfer_gnum(d, _s_elt_gnum);

  s_elt_gnum = NULL;
  _s_elt_gnum = NULL;

  cs_lnum_t num_id = 0;

  do {   /* Loop on sections with equivalent MED element type */

    const fvm_nodal_section_t *section = current_section->section;
    const cs_lnum_t  n_elts_section = section->n_elements;

    if (section->type == current_section->type) {

      /* Ordinary section */

      const cs_lnum_t *vertex_num = section->vertex_num;

      for (elt_id = 0; elt_id < section->n_elements; elt_id++) {
        for (l_id = 0; l_id < stride; l_id++)
          part_vertex_num[num_id++]
            = g_vtx_num[vertex_num[elt_id * stride + vertex_order[l_id]] - 1];
      }

    }
    else {

      /* Tesselated section */

      cs_gnum_t   n_g_sub_elts = 0;
      const cs_lnum_t n_sub_elts
        = fvm_tesselation_n_sub_elements(section->tesselation,
                                         current_section->type);
      const cs_lnum_t *sub_elt_index
        = fvm_tesselation_sub_elt_index(section->tesselation,
                                        current_section->type);
      fvm_tesselation_get_global_size(section->tesselation,
                                      current_section->type,
                                      &n_g_sub_elts,
                                      NULL);

      /* Decode connectivity */

      assert(sub_elt_index[n_elts_section] == n_sub_elts);

      if (n_sub_elts > 0) {

        cs_lnum_t buffer_size = n_sub_elts*stride;
        cs_gnum_t *sub_elt_vtx_gnum = NULL;

        BFT_MALLOC(sub_elt_vtx_gnum, buffer_size, cs_gnum_t);

        fvm_tesselation_decode_g(section->tesselation,
                                 current_section->type,
                                 mesh->global_vertex_num,
                                 current_section->extra_vertex_base,
                                 sub_elt_vtx_gnum);

        /* Convert FVM connectivity to MED connectivity */

        for (elt_id = 0; elt_id < n_sub_elts; elt_id++) {
          for (l_id = 0; l_id < stride; l_id++)
            part_vertex_num[num_id++]
              = sub_elt_vtx_gnum[  (elt_id * stride)
                                 + vertex_order[l_id]];
        }

        BFT_FREE(sub_elt_vtx_gnum);

      }

    } /* End of tesselated section */

    current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous);

  cs_part_to_block_copy_array(d,
                              _med_int_datatype(),
                              stride,
                              part_vertex_num,
                              block_vertex_num);

  cs_part_to_block_destroy(&d);

  BFT_FREE(part_vertex_num);

  /* Write buffers into MED file */
  /*-----------------------------*/

  /* Write connectivity */

  if (w->block_comm != MPI_COMM_NULL) { /* Parallel IO */

    med_int count = (bi.gnum_range[1] > bi.gnum_range[0]) ? 1 : 0;

    med_filter filter = MED_FILTER_INIT;
    retval = MEDfilterBlockOfEntityCr(w->fid,
                                      n_g_elts,
                                      1,
                                      stride,
                                      MED_ALL_CONSTITUENT,
                                      MED_FULL_INTERLACE,
                                      MED_COMPACT_STMODE,
                                      MED_NO_PROFILE,
                                      bi.gnum_range[0], /* start */
                                      n_block_elts,     /* stride */
                                      count,
                                      n_block_elts,     /* blocksize */
                                      0,                /* lastblocksize */
                                      &filter);

    if (retval < 0)
      bft_error(__FILE__, __LINE__, 0,
                _("MEDfilterBlockOfEntityCr() failed for connectivty.\n"
                  "Associated writer: \"%s\"\n"
                  "Associated med_mesh: \"%s\"\n"),
                w->name, med_mesh->name);

    retval = MEDmeshElementConnectivityAdvancedWr(w->fid,
                                                  med_mesh->name,
                                                  MED_NO_DT,
                                                  MED_NO_IT,
                                                  0.0,
                                                  MED_CELL,
                                                  med_section_type,
                                                  MED_NODAL,
                                                  &filter,
                                                  block_vertex_num);

    if (retval < 0)
      bft_error
        (__FILE__, __LINE__, 0,
         _("%s failed to write connectivity:\n"
           "Associated writer: \"%s\"\n"
           "Associated med_mesh_name: \"%s\"\n"
           "Associated MED geometrical element: \"%i\"\n"),
         "MEDmeshElementConnectivityAdvancedWr",
         w->name, med_mesh->name, med_section_type);

    MEDfilterClose(&filter);

  }
  else if (w->rank == 0) { /* Serial IO */

    retval = MEDmeshElementConnectivityWr(w->fid,
                                          med_mesh->name,
                                          MED_NO_DT,
                                          MED_NO_IT,
                                          0.0,
                                          MED_CELL,
                                          med_section_type,
                                          MED_NODAL,
                                          MED_FULL_INTERLACE,
                                          (med_int)n_block_elts,
                                          block_vertex_num);

    if (retval < 0)
      bft_error
        (__FILE__, __LINE__, 0,
         _("%s failed to write connectivity:\n"
           "Associated writer: \"%s\"\n"
           "Associated med_mesh_name: \"%s\"\n"
           "Associated MED geometrical element: \"%i\"\n"),
         "MEDmeshElementConnectivityWr",
         w->name, med_mesh->name, med_section_type);

  } /* If rank == 0 */

  BFT_FREE(block_vertex_num);

  /* Write family numbers */

  _export_families_g(export_sections, w, med_mesh);

  return current_section;
}

#endif /* HAVE_MPI */

/*----------------------------------------------------------------------------
 * Write element family values of a nodal mesh to a MED file in serial mode.
 *
 * parameters:
 *   export_section <-- pointer to MED section helper structure
 *   writer         <-- pointer to associated writer.
 *   med_mesh       <-- pointer to MED mesh structure.
 *
 * returns:
 *  pointer to next EnSight section helper structure in list
 *----------------------------------------------------------------------------*/

static const fvm_writer_section_t *
_export_families_l(const fvm_writer_section_t  *export_section,
                   fvm_to_med_writer_t         *writer,
                   fvm_to_med_mesh_t           *med_mesh)
{
  cs_lnum_t   i;
  med_geometry_type  med_section_type;

  cs_lnum_t   start_id = 0;
  cs_lnum_t   n_elements = 0;

  med_err     retval = 0;
  med_int    *elt_family = NULL;

  const fvm_writer_section_t  *current_section = NULL;

  cs_datatype_t  gc_id_type = CS_DATATYPE_NULL;
  cs_datatype_t  med_family_type = CS_DATATYPE_NULL;

  /* Initialize datatypes */

  switch(sizeof(int)) {
  case 4:
    gc_id_type = CS_INT32;
    break;
  case 8:
    gc_id_type = CS_INT64;
    break;
  default:
    assert(0);
  }

  switch(sizeof(med_int)) {
  case 4:
    med_family_type = CS_INT32;
    break;
  case 8:
    med_family_type = CS_INT64;
    break;
  default:
    assert(0);
  }

  /* Loop on sections to count output size */

  current_section = export_section;
  do {

    const fvm_nodal_section_t  *section = current_section->section;

    if (current_section->type == section->type)
      n_elements += section->n_elements;
    else
      n_elements += fvm_tesselation_n_sub_elements(section->tesselation,
                                                   current_section->type);

    current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous == true);

  BFT_MALLOC(elt_family, n_elements, med_int);

  /* Distribute partition to block values */

  /* loop on sections which should be appended */

  start_id = 0;

  current_section = export_section;
  do {

    const fvm_nodal_section_t  *section = current_section->section;
    const void *src_data[1] = {section->gc_id};

    if (section->gc_id != NULL) {

      fvm_convert_array(1, /* src_dim */
                        0, /* src_dim_shift */
                        1, /* dest_dim */
                        0,
                        section->n_elements,
                        CS_NO_INTERLACE,
                        gc_id_type,
                        med_family_type,
                        0, /* n_parent_lists */
                        NULL,
                        NULL,
                        src_data,
                        elt_family + start_id);

      for (i = 0; i < section->n_elements; i++)
        elt_family[start_id + i] *= -1;

    }
    else {

      for (i = 0; i < section->n_elements; i++)
        elt_family[start_id + i] = 0;

    }

    /* Duplicate values in case of tesselation */

    if (current_section->type != section->type)
      fvm_tesselation_distribute(section->tesselation,
                                 current_section->type,
                                 0,
                                 section->n_elements,
                                 sizeof(med_int),
                                 elt_family + start_id);

    if (current_section->type == section->type)
      start_id += section->n_elements;
    else
      start_id += fvm_tesselation_n_sub_elements(section->tesselation,
                                                 current_section->type);

    current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous == true);

  /* Write values */

  med_section_type = _get_med_elt_type(export_section->type);

  retval = MEDmeshEntityFamilyNumberWr(writer->fid,
                                       med_mesh->name,
                                       MED_NO_DT,
                                       MED_NO_IT,
                                       MED_CELL,
                                       med_section_type,
                                       n_elements,
                                       elt_family);

  if (retval < 0)
    bft_error
      (__FILE__, __LINE__, 0,
       _("MEDmeshEntityFamilyNumberWr() failed to write family numbers:\n"
         "Associated writer: \"%s\"\n"
         "Associated med_mesh_name: \"%s\"\n"
         "Associated MED geometrical element: \"%i\"\n"),
       writer->name, med_mesh->name, med_section_type);

  BFT_FREE(elt_family);

  return current_section;
}

/*----------------------------------------------------------------------------
 * Write strided elements connectivity to a MED file in serial mode
 *
 * parameters:
 *   export_sections       <-- pointer to sections list to export.
 *   writer                <-- pointer to associated writer.
 *   mesh                  <-- pointer to FVM mesh structure.
 *   med_mesh              <-- pointer to MED mesh structure.
 *   export_connect        <-- buffer to export connectivity.
 *
 * returns:
 *  pointer to next MED section structure in list
 *----------------------------------------------------------------------------*/

static const fvm_writer_section_t *
_export_connect_l(const fvm_writer_section_t  *export_sections,
                  fvm_to_med_writer_t         *writer,
                  fvm_to_med_mesh_t           *med_mesh,
                  med_int                     *med_export_connect)
{
  int vertex_order[FVM_MED_MAX_N_NODES];
  cs_lnum_t   i_vtx, i_elt, i_num;
  med_geometry_type  med_section_type;

  med_int  n_export_elements = 0;

  const int stride = fvm_nodal_n_vertices_element[export_sections->type];
  const fvm_writer_section_t *current_section = NULL;

  med_err  retval = 0;

  current_section = export_sections;

  med_section_type = _get_med_elt_type(current_section->type);

  /* Get MED element type vertex order */

  _get_vertex_order(med_section_type,
                    vertex_order);

  /* Gather connectivity from sections sharing the same element type */
  /*-----------------------------------------------------------------*/

  do {   /* Loop on sections with equivalent MED element type */

    const fvm_nodal_section_t *section = current_section->section;

    if (section->type == current_section->type) {

      /* Ordinary section */
      /*------------------*/

      const cs_lnum_t *vertex_num = section->vertex_num;

      /* Convert FVM connectivity to be congruent with MED standard */

      i_num = n_export_elements * stride;
      for (i_elt = 0; i_elt < section->n_elements; i_elt++) {
        for (i_vtx = 0; i_vtx < stride; i_vtx++)
          med_export_connect[i_num++] =
            (med_int)vertex_num[i_elt * stride + vertex_order[i_vtx]];
      }

      n_export_elements +=  section->n_elements;

    }
    else {

      /* Tesselated section */
      /*--------------------*/

      size_t  buffer_size = 0;
      cs_lnum_t   start_id = 0, end_id = 0;
      cs_lnum_t   n_sub_elts = 0, n_sub_loc = 0, n_sub_elts_max = 0;

      cs_lnum_t   *sub_elt_vertex_num = NULL;

      const cs_lnum_t *sub_elt_index
        = fvm_tesselation_sub_elt_index(section->tesselation,
                                        current_section->type);

      n_sub_elts = fvm_tesselation_n_sub_elements(section->tesselation,
                                                  current_section->type);

      fvm_tesselation_get_global_size(section->tesselation,
                                      current_section->type,
                                      NULL,
                                      &n_sub_elts_max);

      buffer_size = CS_MAX(n_sub_elts_max, n_sub_elts/4 + 1);
      buffer_size = CS_MAX(256, buffer_size);

      BFT_MALLOC(sub_elt_vertex_num, buffer_size * stride, cs_lnum_t);

      do {

        end_id
          = fvm_tesselation_decode(section->tesselation,
                                   current_section->type,
                                   start_id,
                                   buffer_size,
                                   current_section->extra_vertex_base,
                                   sub_elt_vertex_num);

        /* Convert FVM connectivity to be congruent with MED standard */

        n_sub_loc = sub_elt_index[end_id] - sub_elt_index[start_id];

        i_num = n_export_elements * stride;
        for (i_elt = 0; i_elt < n_sub_loc; i_elt++) {
          for (i_vtx = 0; i_vtx < stride; i_vtx++)
            med_export_connect[i_num++]
              = (med_int)sub_elt_vertex_num[  (i_elt * stride)
                                            + vertex_order[i_vtx]];
        }

        n_export_elements += n_sub_loc;
        start_id = end_id;

      } while (end_id < section->n_elements);

      BFT_FREE(sub_elt_vertex_num);

    } /* End of tesselated section */

    current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous);

  /* Write buffers into MED file */
  /*-----------------------------*/

  retval = MEDmeshElementConnectivityWr(writer->fid,
                                        med_mesh->name,
                                        MED_NO_DT,
                                        MED_NO_IT,
                                        0.0,
                                        MED_CELL,
                                        med_section_type,
                                        MED_NODAL,
                                        MED_FULL_INTERLACE,
                                        n_export_elements,
                                        med_export_connect);

  if (retval < 0)
    bft_error
      (__FILE__, __LINE__, 0,
       _("MEDmeshElementConnectivityWr() failed to write connectivity:\n"
         "Associated writer: \"%s\"\n"
         "Associated med_mesh_name: \"%s\"\n"
         "Associated MED geometrical element: \"%i\"\n"),
       writer->name, med_mesh->name, med_section_type);

  /* Write family numbers */

  _export_families_l(export_sections, writer, med_mesh);

  return current_section;
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Write polygonal elements connectivity to a MED file in parallel mode
 *
 * parameters:
 *   export_sections  <-- pointer to sections list to export.
 *   w                <-- pointer to associated writer.
 *   mesh             <-- pointer to FVM mesh structure.
 *   med_mesh         <-- pointer to MED mesh structure.
 *----------------------------------------------------------------------------*/

static const fvm_writer_section_t *
_export_nodal_polygons_g(const fvm_writer_section_t  *export_sections,
                         fvm_to_med_writer_t         *w,
                         const fvm_nodal_t           *mesh,
                         fvm_to_med_mesh_t           *med_mesh)
{
  const fvm_writer_section_t  *current_section = NULL;
  const cs_gnum_t *g_vtx_num
    = fvm_io_num_get_global_num(mesh->global_vertex_num);

  med_err  retval = 0;

  /* Get MED element type */

  assert(_get_med_elt_type(export_sections->type) == MED_POLYGON);
  assert(w->discard_polygons == false);

  /* Count elements and prepare block distribution */

  cs_lnum_t   n_part_elts = 0;
  cs_gnum_t   n_g_elts = 0;

  _count_connect_g(export_sections, &n_part_elts, &n_g_elts);

  /* No parallel IO yet for polygons, so use local options */

  const cs_block_dist_info_t bi = cs_block_dist_compute_sizes(w->rank,
                                                              w->n_ranks,
                                                              w->n_ranks,
                                                              w->min_block_size,
                                                              n_g_elts);

  const cs_lnum_t n_block_elts = bi.gnum_range[1] - bi.gnum_range[0];

  const cs_gnum_t *s_elt_gnum
    = fvm_io_num_get_global_num(export_sections->section->global_element_num);
  cs_gnum_t *_s_elt_gnum = _section_elt_gnum(export_sections);

  if (_s_elt_gnum != NULL)
    s_elt_gnum = _s_elt_gnum;

  cs_part_to_block_t *d
    = cs_part_to_block_create_by_gnum(w->comm, bi, n_part_elts, s_elt_gnum);

  if (_s_elt_gnum != NULL)
    cs_part_to_block_transfer_gnum(d, _s_elt_gnum);

  s_elt_gnum = NULL;
  _s_elt_gnum = NULL;

  /* Build global polygon -> vertices index */
  /*----------------------------------------*/

  current_section = export_sections;

  cs_lnum_t *block_index = NULL;
  cs_lnum_t *_part_index = NULL;
  const cs_lnum_t *part_index = current_section->section->vertex_index;

  BFT_MALLOC(block_index, n_block_elts + 1, cs_lnum_t);

  /* Build copy if multiple sections need to be appended,
     point to index otherwise */

  if (n_part_elts > current_section->section->n_elements) {

    cs_lnum_t num_id = 0;

    current_section = export_sections;

    BFT_MALLOC(_part_index, n_part_elts + 1, cs_lnum_t);
    part_index = _part_index;

    _part_index[0] = 0;

    do {   /* Loop on sections with equivalent MED element type */

      const fvm_nodal_section_t *section = current_section->section;
      const cs_lnum_t  n_elts_section = section->n_elements;
      const cs_lnum_t *vertex_index = section->vertex_index;

      for (cs_lnum_t elt_id = 0; elt_id < n_elts_section; elt_id++) {
        _part_index[num_id+1] =   part_index[num_id]
                                + vertex_index[elt_id+1] - vertex_index[elt_id];
        num_id += 1;
      }

      current_section = current_section->next;

    } while (   current_section != NULL
             && current_section->continues_previous);

  }

  /* Build as count, will convert to index later */

  cs_part_to_block_copy_index(d,
                              part_index,
                              block_index);

  /* Gather connectivity from sections sharing the same element type */
  /*-----------------------------------------------------------------*/

  current_section = export_sections;

  cs_lnum_t section_shift = 0;

  med_int *part_connect, *block_connect;

  cs_lnum_t block_size = block_index[n_block_elts];

  BFT_MALLOC(block_connect, block_size, med_int);
  BFT_MALLOC(part_connect, part_index[n_part_elts], med_int);

  do {   /* Loop on sections with equivalent MED element type */

    const fvm_nodal_section_t *section = current_section->section;
    const cs_lnum_t  section_size = section->vertex_index[section->n_elements];

    for (cs_lnum_t v_id = 0; v_id < section_size; v_id++)
      part_connect[v_id + section_shift]
        = g_vtx_num[section->vertex_num[v_id] - 1];

    section_shift += section_size;

    current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous);

  cs_part_to_block_copy_indexed(d,
                                _med_int_datatype(),
                                part_index,
                                part_connect,
                                block_index,
                                block_connect);

  cs_part_to_block_destroy(&d);

  BFT_FREE(part_connect);
  BFT_FREE(_part_index);

  /* Build global block index */

  cs_gnum_t block_end = 0, _block_size = block_size;

  MPI_Scan(&_block_size, &block_end, 1, CS_MPI_GNUM, MPI_SUM, w->comm);
  block_end += 1;
  const cs_gnum_t block_start = block_end - block_size;

  med_int *g_block_index;
  BFT_MALLOC(g_block_index, n_block_elts + 1, med_int);

  for (cs_lnum_t v_id = 0; v_id < n_block_elts+1; v_id++)
    g_block_index[v_id] = block_index[v_id] + block_start;

  BFT_FREE(block_index);

  /* Write buffers into MED file */
  /*-----------------------------*/

  if (w->block_comm != MPI_COMM_NULL) {
    _med_file_close(w);
    _med_file_open_serial(w);
  }

  if (w->rank == 0) {

    /* Write polygonal connectivity into MED file */

    retval = MEDmeshPolygonWr(w->fid,
                              med_mesh->name,
                              MED_NO_DT,
                              MED_NO_IT,
                              0.0,
                              MED_CELL,
                              MED_NODAL,
                              n_block_elts + 1,
                              g_block_index,
                              block_connect);
    if (retval < 0)
      bft_error(__FILE__, __LINE__, 0,
                _("MEDmeshPolygonWr() failed to write connectivity:\n"
                  "Associated writer: \"%s\"\n"
                  "Associated med_mesh_name: \"%s\"\n"),
                w->name, med_mesh->name);

  } /* rank == 0 */

  if (w->block_comm != MPI_COMM_NULL) {
    _med_file_close(w);
    _med_file_open(w, MED_ACC_RDWR);
  }

  BFT_FREE(g_block_index);
  BFT_FREE(block_connect);

  return current_section;
}

#endif /* HAVE_MPI */

/*----------------------------------------------------------------------------
 * Write polygonal element connectivity to a MED file in serial mode
 *
 * parameters:
 *   export_sections       <-- pointer to sections list to export.
 *   writer                <-- pointer to associated writer.
 *   mesh                  <-- pointer to FVM mesh structure.
 *   med_mesh              <-- pointer to MED mesh structure.
 *   med_export_connect    <-- buffer to export connectivity.
 *----------------------------------------------------------------------------*/

static const fvm_writer_section_t *
_export_nodal_polygons_l(const fvm_writer_section_t  *export_sections,
                         fvm_to_med_writer_t         *writer,
                         fvm_to_med_mesh_t           *med_mesh,
                         med_int                     *med_export_connect)
{
  int i_count;
  cs_gnum_t   i;

  int   n_passes = 0;
  cs_gnum_t   n_export_connect = 0;
  cs_gnum_t   n_export_elements = 0;
  med_int   offset = 1;
  med_int  *med_global_vtx_idx = NULL;

  const fvm_writer_section_t *current_section = NULL;

  med_err  retval = 0;

  current_section = export_sections;

  /* Get MED element type */

  assert(_get_med_elt_type(current_section->type) == MED_POLYGON);
  assert(writer->discard_polygons == false);

  /* Gather connectivity from sections sharing the same element type */
  /*-----------------------------------------------------------------*/

  do {   /* Loop on sections with equivalent MED element type */

    const fvm_nodal_section_t *section = current_section->section;
    const cs_gnum_t n_connect_size = (cs_gnum_t)section->connectivity_size;
    const cs_gnum_t n_elements_section
      = fvm_nodal_section_n_g_elements(section);

      n_passes++;

      /* Allocate global buffers */

      BFT_REALLOC(med_global_vtx_idx,
                  n_export_elements + n_elements_section + 1,
                  med_int);

      /*
        Build global vertex connectivity. First, we have to create a global
        index: global_vertex_index. Then, we gather local vertex connectivity
        into med_export_connect. We have to use an offset of n_export_elements
        as there can be several sections of the same type used in fvm_nodal_t.
      */

      for (i = 0; i < n_elements_section + 1; i++)
        med_global_vtx_idx[i + n_export_elements]
          = (med_int)section->vertex_index[i];

      for (i = 0; i < n_connect_size; i++)
        med_export_connect[i + n_export_connect]
          = (med_int)section->vertex_num[i];

      n_export_connect += n_connect_size;
      n_export_elements += n_elements_section + 1;

      current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous);

  assert(n_passes >= 1);

  /* Concatenate med_global_vtx_idx */

  i_count = 0;
  n_export_elements = n_export_elements + 1 - n_passes;
  med_global_vtx_idx[0] = 1;

  for (i = 1; i < n_export_elements; i++) {

    if (med_global_vtx_idx[i + i_count] == 0) {
      i_count++;
      offset = med_global_vtx_idx[i - 1];
    }
    med_global_vtx_idx[i] = med_global_vtx_idx[i + i_count] + offset;
  }

  assert(n_passes == i_count + 1);

  /* Write buffers into MED file */
  /*-----------------------------*/

  /* Write polygonal connectivity into MED file */

  retval = MEDmeshPolygonWr(writer->fid,
                            med_mesh->name,
                            MED_NO_DT,
                            MED_NO_IT,
                            0.0,
                            MED_CELL,
                            MED_NODAL,
                            (med_int)n_export_elements,
                            med_global_vtx_idx,
                            med_export_connect);
  if (retval < 0)
    bft_error(__FILE__, __LINE__, 0,
              _("MEDmeshPolygonWr() failed to write connectivity:\n"
                "Associated writer: \"%s\"\n"
                "Associated med_mesh_name: \"%s\"\n"),
              writer->name, med_mesh->name);

  BFT_FREE(med_global_vtx_idx);

  /* Write family numbers */

  _export_families_l(export_sections, writer, med_mesh);

  return current_section;
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Write polyhedral elements connectivity to a MED file in parallel mode
 *
 * parameters:
 *   export_sections       <-- pointer to sections list to export.
 *   w                     <-- pointer to associated writer.
 *   mesh                  <-- pointer to FVM mesh structure.
 *   med_mesh              <-- pointer to MED mesh structure.
 *----------------------------------------------------------------------------*/

static const fvm_writer_section_t *
_export_nodal_polyhedra_g(const fvm_writer_section_t  *export_sections,
                          fvm_to_med_writer_t         *w,
                          const fvm_nodal_t           *mesh,
                          fvm_to_med_mesh_t           *med_mesh)
{
  const cs_gnum_t *g_vtx_num
    = fvm_io_num_get_global_num(mesh->global_vertex_num);

  med_err  retval = 0;

  /* Get MED element type */

  assert(_get_med_elt_type(export_sections->type) == MED_POLYHEDRON);
  assert(w->discard_polyhedra == false);

  /* Count elements and prepare block distribution */

  cs_lnum_t   n_part_elts = 0;
  cs_gnum_t   n_g_elts = 0;

  _count_connect_g(export_sections, &n_part_elts, &n_g_elts);

  /* No parallel IO yet for polyhedra, so use local options */

  const cs_block_dist_info_t bi = cs_block_dist_compute_sizes(w->rank,
                                                              w->n_ranks,
                                                              w->n_ranks,
                                                              w->min_block_size,
                                                              n_g_elts);

  const cs_lnum_t n_block_elts = bi.gnum_range[1] - bi.gnum_range[0];

  const cs_gnum_t *s_elt_gnum
    = fvm_io_num_get_global_num(export_sections->section->global_element_num);
  cs_gnum_t *_s_elt_gnum = _section_elt_gnum(export_sections);

  if (_s_elt_gnum != NULL)
    s_elt_gnum = _s_elt_gnum;

  cs_part_to_block_t *d
    = cs_part_to_block_create_by_gnum(w->comm, bi, n_part_elts, s_elt_gnum);

  if (_s_elt_gnum != NULL)
    cs_part_to_block_transfer_gnum(d, _s_elt_gnum);

  s_elt_gnum = NULL;
  _s_elt_gnum = NULL;

  /* Build global polyhedron -> faces index */
  /*----------------------------------------*/

  cs_lnum_t *block_f_index = NULL, *_part_f_index = NULL;
  const fvm_writer_section_t  *current_section = export_sections;
  const cs_lnum_t *part_f_index = current_section->section->face_index;

  BFT_MALLOC(block_f_index, n_block_elts + 1, cs_lnum_t);

  /* Build copy if multiple sections need to be appended,
     point to index otherwise */

  if (n_part_elts > current_section->section->n_elements) {

    cs_lnum_t num_id = 0;

    current_section = export_sections;

    BFT_MALLOC(_part_f_index, n_part_elts + 1, cs_lnum_t);
    part_f_index = _part_f_index;

    _part_f_index[0] = 0;

    do {   /* Loop on sections with equivalent MED element type */

      const fvm_nodal_section_t *section = current_section->section;
      const cs_lnum_t *face_index = section->face_index;

      for (cs_lnum_t elt_id = 0; elt_id < section->n_elements; elt_id++) {
        _part_f_index[num_id+1] =   part_f_index[num_id]
                                  + face_index[elt_id+1] - face_index[elt_id];
        num_id += 1;
      }

      current_section = current_section->next;

    } while (   current_section != NULL
             && current_section->continues_previous);

  }

  cs_part_to_block_copy_index(d,
                              part_f_index,
                              block_f_index);

  cs_lnum_t part_f_count = part_f_index[n_part_elts];
  cs_lnum_t block_f_count = block_f_index[n_block_elts];

  /* Build global polyhedron -> vertices index */
  /*-------------------------------------------*/

  cs_lnum_t *block_v_index = NULL, *part_v_index = NULL;
  med_int *block_f_v_index = NULL, *part_f_v_count = NULL;

  BFT_MALLOC(block_v_index, n_block_elts + 1, cs_lnum_t);
  BFT_MALLOC(part_v_index, n_part_elts + 1, cs_lnum_t);

  BFT_MALLOC(block_f_v_index, block_f_count+1, med_int);
  BFT_MALLOC(part_f_v_count, part_f_count, med_int);

  current_section = export_sections;

  part_v_index[0] = 0;

  part_f_count = 0; /* reset, will be re-incremented */

  cs_lnum_t elt_count = 0;

  do {   /* Loop on sections with equivalent MED element type */

    const fvm_nodal_section_t *section = current_section->section;
    const cs_lnum_t *face_index = section->face_index;

    for (cs_lnum_t elt_id = 0; elt_id < section->n_elements; elt_id++) {

      cs_lnum_t elt_v_count = 0;

      for (cs_lnum_t i = face_index[elt_id]; i < face_index[elt_id + 1]; i++) {
        cs_lnum_t face_id = CS_ABS(section->face_num[i]) - 1;
        cs_lnum_t n_f_vertices =   section->vertex_index[face_id + 1]
                                 - section->vertex_index[face_id];
        elt_v_count += n_f_vertices;
        part_f_v_count[part_f_count++] = n_f_vertices;
      }

      part_v_index[elt_count+1] = part_v_index[elt_count] + elt_v_count;

      elt_count += 1;

    }

    current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous);

  cs_part_to_block_copy_index(d,
                              part_v_index,
                              block_v_index);

  /* Distribute face -> vertex counts; on the block distribution,
     it will be later converted to an index, so shift it by one place */

  cs_part_to_block_copy_indexed(d,
                                _med_int_datatype(),
                                part_f_index,
                                part_f_v_count,
                                block_f_index,
                                block_f_v_index + 1);

  BFT_FREE(part_f_v_count);
  BFT_FREE(_part_f_index);
  part_f_index = NULL;

  /* Build connectivity from sections sharing the same element type */
  /*----------------------------------------------------------------*/

  cs_lnum_t part_connect_size = 0;

  med_int *block_connect = NULL, *part_connect = NULL;

  BFT_MALLOC(block_connect, block_v_index[n_block_elts], med_int);
  BFT_MALLOC(part_connect, part_v_index[n_part_elts], med_int);

  current_section = export_sections;

  do {   /* Loop on sections with equivalent MED element type */

    const fvm_nodal_section_t *section = current_section->section;
    const cs_lnum_t *face_index = section->face_index;

    for (cs_lnum_t elt_id = 0; elt_id < section->n_elements; elt_id++) {
      for (cs_lnum_t i = face_index[elt_id]; i < face_index[elt_id + 1]; i++) {

        if (section->face_num[i] > 0) {
          cs_lnum_t face_id = section->face_num[i] - 1;
          for (cs_lnum_t v_id = section->vertex_index[face_id];
               v_id < section->vertex_index[face_id + 1];
               v_id++)
            part_connect[part_connect_size++]
              = g_vtx_num[section->vertex_num[v_id] - 1];
        }
        else { /* face_num < 0 */
          cs_lnum_t face_id = -section->face_num[i] - 1;
          for (cs_lnum_t v_id = section->vertex_index[face_id + 1] - 1;
               v_id >= section->vertex_index[face_id];
               v_id--)
            part_connect[part_connect_size++]
              = g_vtx_num[section->vertex_num[v_id] - 1];
        }

      }
    }

    current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous);

  cs_part_to_block_copy_indexed(d,
                                _med_int_datatype(),
                                part_v_index,
                                part_connect,
                                block_v_index,
                                block_connect);

  BFT_FREE(part_connect);
  BFT_FREE(part_v_index);

  cs_part_to_block_destroy(&d); /* All distribution done */

  /* Build global block indexes: 0 for faces, 1 for vertices */

  cs_lnum_t block_f_size = block_f_index[n_block_elts];
  cs_lnum_t block_v_size = block_v_index[n_block_elts];

  cs_gnum_t block_end[2] = {0, 0};
  cs_gnum_t block_size[2] = {block_f_size, block_v_size};

  MPI_Scan(&block_size, &block_end, 2, CS_MPI_GNUM, MPI_SUM, w->comm);
  block_end[0] += 1, block_end[1] += 1;

  const cs_gnum_t block_f_start = block_end[0] - block_size[0];
  const cs_gnum_t block_v_start = block_end[1] - block_size[1];

  BFT_FREE(block_v_index);

  /* Convert face vertex counts to index */

  block_f_v_index[0] = block_v_start;
  for (cs_lnum_t i = 0; i < block_f_size; i++)
    block_f_v_index[i+1] += block_f_v_index[i];

  /* Build global cell -> faces index */

  med_int *g_block_f_index;
  BFT_MALLOC(g_block_f_index, n_block_elts + 1, med_int);

  for (cs_lnum_t i = 0; i < n_block_elts+1; i++)
    g_block_f_index[i] = block_f_index[i] + block_f_start;

  BFT_FREE(block_f_index);

  /* Write buffers into MED file */
  /*-----------------------------*/

  if (w->block_comm != MPI_COMM_NULL) {
    _med_file_close(w);
    _med_file_open_serial(w);
  }

  if (w->rank == 0) {

    /* Write polyhedral connectivity into MED file */

    retval = MEDmeshPolyhedronWr(w->fid,
                                 med_mesh->name,
                                 MED_NO_DT,
                                 MED_NO_IT,
                                 0.0,
                                 MED_CELL,
                                 MED_NODAL,
                                 n_block_elts + 1,
                                 g_block_f_index,
                                 block_f_size + 1,
                                 block_f_v_index,
                                 block_connect);

    if (retval < 0)
      bft_error(__FILE__, __LINE__, 0,
                _("MEDmeshPolyhedronWr() failed to write connectivity:\n"
                  "Associated writer: \"%s\"\n"
                  "Associated med_mesh_name: \"%s\"\n"),
                w->name, med_mesh->name);

  } /* rank == 0 */

  if (w->block_comm != MPI_COMM_NULL) {
    _med_file_close(w);
    _med_file_open(w, MED_ACC_RDWR);
  }

  BFT_FREE(g_block_f_index);
  BFT_FREE(block_connect);
  BFT_FREE(block_f_v_index);

  /* Write family numbers */

  _export_families_g(export_sections, w, med_mesh);

  return current_section;
}

#endif /* HAVE_MPI */

/*----------------------------------------------------------------------------
 * Write polyhedral elements connectivity to a MED file in serial mode
 *
 * parameters:
 *   export_sections       <-- pointer to sections list to export.
 *   writer                <-- pointer to associated writer.
 *   mesh                  <-- pointer to FVM mesh structure.
 *   med_mesh              <-- pointer to MED mesh structure.
 *   med_elt_type          <-- MED geometrical element type to export.
 *   n_g_element_section   <-- global number of elements per section.
 *   export_connect        <-- buffer to export connectivity.
 *----------------------------------------------------------------------------*/

static const fvm_writer_section_t *
_export_nodal_polyhedra_l(const fvm_writer_section_t  *export_sections,
                          fvm_to_med_writer_t         *writer,
                          fvm_to_med_mesh_t           *med_mesh,
                          med_int                     *med_export_connect)
{
  cs_lnum_t   i_face_idx, face_id, vtx_id, i_elt;
  cs_gnum_t   i;

  cs_gnum_t   i_cell = 1, i_face = 1, i_connect = 0;
  cs_gnum_t   n_export_elements = 0;
  cs_gnum_t   n_export_faces = 0;

  med_int  *med_face_lengths = NULL;
  med_int  *med_cell_lengths = NULL;

  const fvm_writer_section_t *current_section = NULL;

  med_err  retval = 0;

  current_section = export_sections;

  /* Get MED element type */

  assert(_get_med_elt_type(current_section->type) == MED_POLYHEDRON);
  assert(writer->discard_polyhedra == false);

  /* Gather connectivity from sections sharing the same element type */
  /*-----------------------------------------------------------------*/

  do {   /* Loop on sections with equivalent MED element type */

    const fvm_nodal_section_t *section = current_section->section;

    cs_gnum_t   n_faces_section = 0;

    n_faces_section = section->face_index[section->n_elements];

    BFT_REALLOC(med_face_lengths,
                n_faces_section + 1 + n_export_faces,
                med_int);

    BFT_REALLOC(med_cell_lengths,
                section->n_elements + 1 + n_export_elements,
                med_int);

    /*
      Build locally:
      - med_face_lengths (number of vertices per face),
      - med_cell_lengths (number of faces per cell) and
      - cells -> vertices connectivity for each polyhedral section
    */

    for (i_elt = 0; i_elt < section->n_elements; i_elt++) {

      med_cell_lengths[i_cell++] = (  (med_int)section->face_index[i_elt + 1]
                                    - (med_int)section->face_index[i_elt]);

      for (i_face_idx = section->face_index[i_elt];
           i_face_idx < section->face_index[i_elt + 1];
           i_face_idx++) {

        face_id = CS_ABS(section->face_num[i_face_idx]) - 1;

        med_face_lengths[i_face] = (  (med_int)section->vertex_index[face_id + 1]
                                    - (med_int)section->vertex_index[face_id]);
        i_face++;

        if (section->face_num[i_face_idx] > 0) {
          for (vtx_id = section->vertex_index[face_id];
               vtx_id < section->vertex_index[face_id + 1];
               vtx_id++)
            med_export_connect[i_connect++] =
              (med_int)section->vertex_num[vtx_id];
          }
          else { /* face_num < 0 */
            for (vtx_id = section->vertex_index[face_id + 1] - 1;
                 vtx_id >= section->vertex_index[face_id];
                 vtx_id--)
              med_export_connect[i_connect++] =
                (med_int)section->vertex_num[vtx_id];
          }

      } /* End of loop on faces */

    } /* End of loop on cells */

    n_export_elements += section->n_elements;
    n_export_faces += n_faces_section;

    current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous);

  /* Write buffers into MED file */
  /*-----------------------------*/

  /* Create MED cells -> faces index from cell_lengths. */

  med_cell_lengths[0] = 1;

  for (i = 1; i < n_export_elements + 1; i++)
    med_cell_lengths[i] = med_cell_lengths[i] + med_cell_lengths[i-1];

  /* Create MED faces -> vertices index to export from face_lengths */

  med_face_lengths[0] = 1;

  for (i = 1; i < n_export_faces + 1; i++)
    med_face_lengths[i] = med_face_lengths[i] + med_face_lengths[i-1];

  /* Write polyhedral connectivity into MED file */

  retval = MEDmeshPolyhedronWr(writer->fid,
                               med_mesh->name,
                               MED_NO_DT,
                               MED_NO_IT,
                               0.0,
                               MED_CELL,
                               MED_NODAL,
                               (med_int)n_export_elements + 1,
                               med_cell_lengths,
                               (med_int)n_export_faces + 1,
                               med_face_lengths,
                               med_export_connect);

  if (retval < 0)
    bft_error(__FILE__, __LINE__, 0,
              _("MEDmeshPolyhedronWr() failed to write connectivity:\n"
                "Associated writer: \"%s\"\n"
                "Associated med_mesh_name: \"%s\"\n"),
              writer->name, med_mesh->name);

  BFT_FREE(med_cell_lengths);
  BFT_FREE(med_face_lengths);

  /* Write family numbers */

  _export_families_l(export_sections, writer, med_mesh);

  return current_section;
}

/*----------------------------------------------------------------------------
 * Write "trivial" point elements
 *
 * parameters:
 *   writer          <-- pointer to associated writer.
 *   mesh            <-- pointer to FVM mesh structure.
 *   med_mesh        <-- pointer to MED mesh structure.
 *----------------------------------------------------------------------------*/

static void
_export_point_elements(fvm_to_med_writer_t   *writer,
                       const fvm_nodal_t     *mesh,
                       fvm_to_med_mesh_t     *med_mesh)
{
  cs_gnum_t  i_elt;

  cs_gnum_t  n_g_vertices = fvm_nodal_get_n_g_vertices(mesh);

  med_int   *export_connect = NULL;

  med_err  retval = 0;

  /* TODO add MPI IO version */

  if (writer->rank == 0) {

    /* Prepare buffer */

    BFT_MALLOC(export_connect, n_g_vertices, med_int);

    for (i_elt = 0; i_elt < n_g_vertices; i_elt++)
      export_connect[i_elt] = i_elt + 1;

    /* Write into MED file */
    /*---------------------*/

    retval = MEDmeshElementConnectivityWr(writer->fid,
                                          med_mesh->name,
                                          MED_NO_DT,
                                          MED_NO_IT,
                                          0.0,
                                          MED_CELL,
                                          MED_POINT1,
                                          MED_NODAL,
                                          MED_FULL_INTERLACE,
                                          (med_int)n_g_vertices,
                                          export_connect);

    if (retval < 0)
      bft_error
        (__FILE__, __LINE__, 0,
         _("MEDmeshElementConnectivityWr() failed to write connectivity:\n"
           "Associated writer: \"%s\"\n"
           "Associated med_mesh_name: \"%s\"\n"
           "Associated MED geometrical element: \"%i\"\n"),
         writer->name, med_mesh->name, (int)MED_POINT1);

    BFT_FREE(export_connect);

  } /* If rank == 0 */
}

/*----------------------------------------------------------------------------
 * Output function for field values.
 *
 * This function is passed to fvm_writer_field_helper_output_* functions.
 *
 * parameters:
 *   context      <-> pointer to writer and field context
 *   datatype     <-- output datatype
 *   dimension    <-- output field dimension
 *   component_id <-- output component id (if non-interleaved)
 *   block_start  <-- start global number of element for current block
 *   block_end    <-- past-the-end global number of element for current block
 *   buffer       <-> associated output buffer
 *----------------------------------------------------------------------------*/

static void
_field_output(void           *context,
              cs_datatype_t   datatype,
              int             dimension,
              int             component_id,
              cs_gnum_t       block_start,
              cs_gnum_t       block_end,
              void           *buffer)
{
  CS_UNUSED(datatype);
#if !defined(HAVE_MED_MPI)
  CS_UNUSED(dimension);
#endif
  CS_UNUSED(component_id);

  bool serial_io = true;

  med_err retval = 0;

  _med_context_t *c = context;

  fvm_to_med_writer_t  *w = c->writer;

  med_int block_sub_size = block_end - block_start;

  cs_gnum_t n_g_elements = c->n_g_elts;

  if (block_sub_size < 0)
    block_sub_size = 0;

#if defined(HAVE_MED_MPI)

  if (w->block_comm != MPI_COMM_NULL) { /* Parallel IO */

    serial_io = false;

    med_int count = (block_sub_size > 0) ? 1 : 0;

    med_filter filter = MED_FILTER_INIT;
    retval = MEDfilterBlockOfEntityCr(w->fid,
                                      n_g_elements,
                                      1,
                                      dimension,
                                      MED_ALL_CONSTITUENT,
                                      MED_FULL_INTERLACE,
                                      MED_COMPACT_STMODE,
                                      MED_NO_PROFILE,
                                      block_start,     /*start */
                                      block_sub_size,  /*stride */
                                      count,
                                      block_sub_size,  /* blocksize */
                                      0,               /* lastblocksize */
                                      &filter);

    if (retval < 0)
      bft_error(__FILE__, __LINE__, 0,
                _("MEDfilterBlockOfEntityCr() failed for field values.\n"
                  "Associated writer: \"%s\"\n"
                  "Associated med_mesh: \"%s\"\n"),
                w->name, c->mesh_name);

    retval = MEDfieldValueAdvancedWr(w->fid,
                                     c->field_name,
                                     c->time_step,
                                     MED_NO_IT,
                                     c->time_value,
                                     c->entity_type,
                                     c->section_type,
                                     MED_NO_LOCALIZATION,
                                     &filter,
                                     (const unsigned char *)buffer);

    MEDfilterClose(&filter);

  }

#endif /* defined(HAVE_MED_MPI) */

  if (w->rank == 0 && serial_io) {  /* Output using a single rank */

    retval = MEDfieldValueWr(w->fid,
                             c->field_name,
                             c->time_step,
                             MED_NO_IT,
                             c->time_value,
                             c->entity_type,
                             c->section_type,
                             MED_FULL_INTERLACE,
                             MED_ALL_CONSTITUENT,
                             block_sub_size,
                             (const unsigned char *)buffer);

  }

  if (retval < 0)
    bft_error(__FILE__, __LINE__, 0,
              "_field_output() failed to write field values\n"
              "Associated fieldname: \"%s\"\n"
              "Associated med mesh: \"%s\"\n"
              "Associated writer name: \"%s\"\n",
              c->field_name, c->mesh_name, w->name);
}

/*----------------------------------------------------------------------------
 * Write field values associated with element values of a nodal mesh to
 * a MED file.
 *
 * Output fields are interlaced. Input arrays may be interlaced or not.
 *
 * parameters:
 *   export_section   <-- pointer to MED section helper structure
 *   w                <-- pointer to writer structure
 *   helper           <-- pointer to general writer helper structure
 *   dim              <-- field dimension
 *   interlace        <-- indicates if field in memory is interlaced
 *   n_parent_lists   <-- indicates if field values are to be obtained
 *                        directly through the local entity index (when 0) or
 *                        through the parent entity numbers (when 1 or more)
 *   parent_num_shift <-- parent list to common number index shifts;
 *                        size: n_parent_lists
 *   datatype         <-- indicates the data type of (source) field values
 *   field_values     <-- array of associated field value arrays
 *   med_mesh_name    <-- MED name of the mesh on which the field is defined
 *   med_field_name   <-- MED name of the field to export.
 *   time_step        <-- number of the current time step
 *   time_value       <-- associated time value
 *
 * returns:
 *  pointer to next section helper structure in list
 *----------------------------------------------------------------------------*/

static const fvm_writer_section_t *
_export_field_values_e(const fvm_writer_section_t      *export_section,
                       fvm_to_med_writer_t             *w,
                       fvm_writer_field_helper_t       *helper,
                       int                              dim,
                       cs_interlace_t                   interlace,
                       int                              n_parent_lists,
                       const cs_lnum_t                  parent_num_shift[],
                       cs_datatype_t                    datatype,
                       const void                *const field_values[],
                       char                            *med_mesh_name,
                       char                            *med_field_name,
                       int                              time_step,
                       double                           time_value)
{
  _med_context_t c;

  c.writer = w;
  c.mesh_name = med_mesh_name;
  c.field_name = med_field_name;
  c.entity_type = MED_CELL;

  c.section_type = _get_med_elt_type(export_section->type);
  c.time_step = time_step;
  c.time_value = time_value;

  cs_lnum_t   n_part_elts = 0; /* dummy argument */
  _count_connect_g(export_section, &n_part_elts, &(c.n_g_elts));

  return fvm_writer_field_helper_output_e(helper,
                                          &c,
                                          export_section,
                                          dim,
                                          interlace,
                                          NULL, /* component order */
                                          n_parent_lists,
                                          parent_num_shift,
                                          datatype,
                                          field_values,
                                          _field_output);
}

/*----------------------------------------------------------------------------
 * Write per node field values into MED writer.
 *
 * Output fields are non interlaced. Input arrays may be interlaced or not.
 *
 * parameters:
 *   export_section     <-- pointer to section helper structure
 *   w                  <-- pointer to associated writer.
 *   helper             <-- pointer to general writer helper structure
 *   dimension          <-- input field dimension
 *   interlace          <-- indicates if field in memory is interlaced
 *   n_parent_lists     <-- indicates if field values are to be obtained
 *                          directly through the local entity index (when 0) or
 *                          through the parent entity numbers (when 1 or more)
 *   parent_num_shift   <-- parent list to common number index shifts;
 *                          size: n_parent_lists
 *   datatype           <-- indicates the data type of (source) field values
 *   field_values       <-- array of associated field value arrays
 *   med_mesh_name      <-- MED name of the mesh on which the field is defined
 *   med_field_name     <-- MED name of the field to export.
 *   time_step          <-- number of the current time step
 *   time_value         <-- associated time value
 *----------------------------------------------------------------------------*/

static void
_export_field_values_n(const fvm_nodal_t               *mesh,
                       fvm_to_med_writer_t             *w,
                       fvm_writer_field_helper_t       *helper,
                       int                              dimension,
                       cs_interlace_t                   interlace,
                       int                              n_parent_lists,
                       const cs_lnum_t                  parent_num_shift[],
                       cs_datatype_t                    datatype,
                       const void                *const field_values[],
                       char                            *med_mesh_name,
                       char                            *med_field_name,
                       int                              time_step,
                       double                           time_value)
{
  /* Get total number of vertices */

  cs_lnum_t   n_extra_vertices;
  cs_gnum_t   n_g_extra_vertices = 0;
  cs_gnum_t    n_g_vertices = (mesh->global_vertex_num != NULL) ?
      fvm_io_num_get_global_count(mesh->global_vertex_num)
    : (cs_gnum_t)(mesh->n_vertices);

  fvm_writer_count_extra_vertices(mesh,
                                  w->divide_polyhedra,
                                  &n_g_extra_vertices,
                                  &n_extra_vertices);

  n_g_vertices += n_g_extra_vertices;

  _med_context_t c;

  c.writer = w;
  c.mesh_name = med_mesh_name;
  c.field_name = med_field_name;
  c.entity_type = MED_NODE;

  c.section_type = MED_POINT1;
  c.time_step = time_step;
  c.time_value = time_value;

  c.n_g_elts = n_g_vertices;

  fvm_writer_field_helper_output_n(helper,
                                   &c,
                                   mesh,
                                   dimension,
                                   interlace,
                                   NULL, /* component order */
                                   n_parent_lists,
                                   parent_num_shift,
                                   datatype,
                                   field_values,
                                   _field_output);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Returns number of library version strings associated with the MED format.
 *
 * The first associated version string should corresponds to the MED library,
 * The second to the HDF5 library.
 *
 * returns:
 *   number of library version strings associated with the MED format.
 *----------------------------------------------------------------------------*/

int
fvm_to_med_n_version_strings(void)
{
  return 2;
}

/*----------------------------------------------------------------------------
 * Returns a library version string associated with the MED format.
 *
 * The first associated version string should corresponds to the MED library,
 * The second to the HDF5 library.
 *
 * In certain cases, when using dynamic libraries, fvm may be compiled
 * with one library version, and linked with another. If both run-time
 * and compile-time version information is available, this function
 * will return the run-time version string by default.
 *
 * Setting the compile_time flag to 1, the compile-time version string
 * will be returned if this is different from the run-time version.
 * If the version is the same, or only one of the 2 version strings are
 * available, a NULL character string will be returned with this flag set.
 *
 * parameters:
 *   string_index <-- index in format's version string list (0 to n-1)
 *   compile_time <-- 0 by default, 1 if we want the compile-time version
 *                    string, if different from the run-time version.
 *
 * returns:
 *   pointer to constant string containing the library's version.
 *----------------------------------------------------------------------------*/

const char *
fvm_to_med_version_string(int string_index,
                          int compile_time_version)
{
  const char * retval = NULL;

  if (compile_time_version) {

    if (string_index == 0) {
      if (MED_NUM_RELEASE > -1)
        snprintf(_med_version_string_[1], 31, "MED %d.%d.%d",
                 MED_MAJOR_NUM, MED_MINOR_NUM, MED_RELEASE_NUM);
      else
        snprintf(_med_version_string_[1], 31, "MED %d.%d",
                 MED_MAJOR_NUM, MED_MINOR_NUM);
      _med_version_string_[1][31] = '\0';
      retval = _med_version_string_[1];
    }
    else if (string_index == 1) {
      snprintf(_hdf5_version_string_[1], 15, "HDF5 %d.%d.%d",
               H5_VERS_MAJOR, H5_VERS_MINOR, H5_VERS_RELEASE);
      _hdf5_version_string_[1][31] = '\0';
      retval = _hdf5_version_string_[1];
    }

  }

  else {
    if (string_index == 0) {
      med_int   med_major, med_minor, med_release;
      MEDlibraryNumVersion(&med_major, &med_minor, &med_release);
      snprintf(_med_version_string_[0], 31, "MED %d.%d.%d",
               (int)med_major, (int)med_minor, (int)med_release);
      _med_version_string_[0][31] = '\0';
      retval = _med_version_string_[0];
    }
    else if (string_index == 1) {
      med_int  h5_major, h5_minor, h5_release;
      MEDlibraryHdfNumVersion(&h5_major, &h5_minor, &h5_release);
      snprintf(_hdf5_version_string_[0], 15, "HDF5 %d.%d.%d",
               (int)h5_major, (int)h5_minor, (int)h5_release);
      _hdf5_version_string_[0][31] = '\0';
      retval = _hdf5_version_string_[0];
    }
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Initialize FVM to MED file writer.
 *
 * Options are:
 *   discard_polygons    do not output polygons or related values
 *   discard_polyhedra   do not output polyhedra or related values
 *   divide_polygons     tesselate polygons with triangles
 *   divide_polyhedra    tesselate polyhedra with tetrahedra and pyramids
 *                       (adding a vertex near each polyhedron's center)
 *   serial_io           force serial IO even when parallel IO is available
 *   update              open file in update mode
 *
 * parameters:
 *   name           <-- base output case name.
 *   options        <-- whitespace separated, lowercase options list
 *   time_dependecy <-- indicates if and how meshes will change with time
 *   comm           <-- associated MPI communicator.
 *
 * returns:
 *   pointer to opaque MED writer structure.
 *----------------------------------------------------------------------------*/

#if defined(HAVE_MPI)
void *
fvm_to_med_init_writer(const char                   *name,
                       const char                   *path,
                       const char                   *options,
                       const fvm_writer_time_dep_t   time_dependency,
                       const MPI_Comm                comm)
#else
void *
fvm_to_med_init_writer(const char                   *name,
                       const char                   *path,
                       const char                   *options,
                       const fvm_writer_time_dep_t   time_dependency)
#endif
{
  fvm_to_med_writer_t  *writer = NULL;

  int  i;
  int  filename_length, name_length, path_length;

  /* Initialize writer */

  BFT_MALLOC(writer, 1, fvm_to_med_writer_t);

  writer->n_med_meshes = 0;
  writer->n_fields  = 0;
  writer->med_meshes   = NULL;
  writer->fields = NULL;

  writer->n_time_steps   = 0;
  writer->time_steps     = NULL;
  writer->time_values    = NULL;
  writer->time_dependency = time_dependency;

  /* Parallel parameters */

  writer->rank = 0;
  writer->n_ranks = 1;

#if defined(HAVE_MPI)

  writer->comm = comm;
  writer->block_comm = MPI_COMM_NULL;
  {
    int mpi_flag, rank, n_ranks;
    MPI_Initialized(&mpi_flag);

    if (mpi_flag && comm != MPI_COMM_NULL) {
      MPI_Comm_rank(writer->comm, &rank);
      MPI_Comm_size(writer->comm, &n_ranks);
      writer->rank = rank;
      writer->n_ranks = n_ranks;
    }
  }

#endif

  writer->min_rank_step = writer->n_ranks;
  writer->min_block_size = 0;

#if defined(HAVE_MED_MPI)

  /* Assign default IO communicator from general settings
     if compatible with current communicator */

  int min_rank_step = 1;
  MPI_Comm w_block_comm, w_comm;
  cs_file_get_default_comm(&min_rank_step, NULL,
                           &w_block_comm, &w_comm);

  if (min_rank_step < writer->min_rank_step) {
    if (comm == w_comm) {
      writer->min_rank_step = min_rank_step;
      writer->block_comm = w_block_comm;
    }
    else {
      writer->min_rank_step = min_rank_step;
      writer->block_comm = comm;
    }
  }

#endif

  /* Reading options */

  writer->allow_update = false;

  writer->discard_polygons = false;
  writer->discard_polyhedra = false;
  writer->divide_polygons = false;
  writer->divide_polyhedra = false;

  if (options != NULL) {

    int i1, i2, l_opt;
    int l_tot = strlen(options);

    i1 = 0; i2 = 0;
    while (i1 < l_tot) {

      for (i2 = i1 ; i2 < l_tot && options[i2] != ' ' ; i2++);
      l_opt = i2 - i1;

      if (        (l_opt == 16)
               && (strncmp(options + i1, "discard_polygons", l_opt) == 0))
        writer->discard_polygons = true;

      else if (   (l_opt == 17)
               && (strncmp(options + i1, "discard_polyhedra", l_opt) == 0))
        writer->discard_polyhedra = true;

      else if (   (l_opt == 15)
               && (strncmp(options + i1, "divide_polygons", l_opt) == 0))
        writer->divide_polygons = true;

      else if (   (l_opt == 16)
               && (strncmp(options + i1, "divide_polyhedra", l_opt) == 0))
        writer->divide_polyhedra = true;

      else if (   (l_opt == 9)
               && (strncmp(options + i1, "serial_io", l_opt) == 0)) {
        writer->min_rank_step = writer->n_ranks;
#if defined(HAVE_MPI)
        writer->block_comm = MPI_COMM_NULL;
#endif
      }

      else if (   (l_opt == 6)
               && (strncmp(options + i1, "update", l_opt) == 0))
        writer->allow_update = true;

      for (i1 = i2 + 1 ; i1 < l_tot && options[i1] == ' ' ; i1++);

    }

  }

  /* Writer name and file name */

  name_length = strlen(name);
  if (name_length == 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Empty MED filename."));

  BFT_MALLOC(writer->name, name_length + 1, char);
  strcpy(writer->name, name);

  for (i = 0; i < name_length; i++) {
    if (writer->name[i] == ' ' || writer->name[i] == '\t')
      writer->name[i] = '_';
  }

  if (path != NULL)
    path_length = strlen(path);
  else
    path_length = 0;
  filename_length = path_length + name_length + strlen(".med");
  BFT_MALLOC(writer->filename, filename_length + 1, char);

  if (path != NULL)
    strcpy(writer->filename, path);
  else
    writer->filename[0] = '\0';

  strcat(writer->filename, writer->name);
  strcat(writer->filename, ".med");

  writer->filename[filename_length] = '\0';
  writer->name[name_length] = '\0';

  /* Open MED file */

  writer->is_open = false;
  writer->fid = -1;

  _med_file_open(writer, MED_ACC_CREAT);

  return writer;
}

/*----------------------------------------------------------------------------
 * Finalize FVM to MED file writer.
 *
 * parameters:
 *   this_writer_p <-- pointer to opaque MED writer structure.
 *
 * returns:
 *   NULL pointer.
 *----------------------------------------------------------------------------*/

void *
fvm_to_med_finalize_writer(void  *this_writer_p)
{
  int i;

  fvm_to_med_writer_t  *writer
                        = (fvm_to_med_writer_t *)this_writer_p;

  assert(writer != NULL);

  /* Close MED File */
  /*----------------*/

  if (writer->is_open == true)
    _med_file_close(writer);
  else
    assert(writer->fid < 0);

  /* Free structures */
  /*-----------------*/

  BFT_FREE(writer->name);
  BFT_FREE(writer->filename);
  BFT_FREE(writer->time_values);
  BFT_FREE(writer->time_steps);

  /* Free fvm_to_med_mesh_t structure */

  for (i = 0; i < writer->n_med_meshes; i++)
    BFT_FREE(writer->med_meshes[i]);
  BFT_FREE(writer->med_meshes);

  for (i = 0; i < writer->n_fields; i++)
    BFT_FREE(writer->fields[i]);

  BFT_FREE(writer->fields);

  /* Free fvm_to_med_writer_t structure */

  BFT_FREE(writer);

  return NULL;
}

/*----------------------------------------------------------------------------
 * Associate new time step with a MED mesh.
 *
 * parameters:
 *   this_writer <-- pointer to associated writer
 *   time_step   <-- time step number
 *   time_value  <-- time_value number
 *----------------------------------------------------------------------------*/

void
fvm_to_med_set_mesh_time(void          *this_writer,
                         int            time_step,
                         double         time_value)
{
  int n_vals;
  fvm_to_med_writer_t  *writer = (fvm_to_med_writer_t *)this_writer;

  const char time_value_err_string[] =
    N_("The time value associated with time step <%d> equals <%g>,\n"
       "but time value <%g> has already been associated with this time step.\n");

  assert(writer != NULL);

  /* First verification on time step */

  if (time_step < 0) {
    if (writer->time_dependency == FVM_WRITER_FIXED_MESH)
      return;
    else
      bft_error(__FILE__, __LINE__, 0,
                _("The given time step value should be >= 0, and not %d\n"),
                time_step);
  }

  if (   writer->time_steps != NULL
      && writer->time_values != NULL) {

    n_vals = writer->n_time_steps;
    if (time_step < writer->time_steps[n_vals - 1])
      bft_error(__FILE__, __LINE__, 0,
                _("The given time step value should be >= %d, and not %d\n"),
                writer->time_steps[n_vals - 1], time_step);

    /* Verifications on time value */

    else if (time_step == writer->time_steps[n_vals - 1]) {
      if (   time_value < writer->time_values[n_vals - 1] - 1.e-16
          || time_value > writer->time_values[n_vals - 1] + 1.e-16)
        bft_error(__FILE__, __LINE__, 0,
                  _(time_value_err_string), time_step,
                  time_value, writer->time_values[n_vals - 1]);
    }
    else { /* Add a new time step and time value */
      writer->n_time_steps += 1;
      n_vals = writer->n_time_steps;

      BFT_REALLOC(writer->time_values, n_vals, double);
      BFT_REALLOC(writer->time_steps, n_vals, int);

      writer->time_values[n_vals - 1] = time_value;
      writer->time_steps[n_vals - 1] = time_step;
    }
  }
  else { /* Setting of the first time step and time value */
    writer->n_time_steps += 1;
    n_vals = writer->n_time_steps;

    BFT_REALLOC(writer->time_values, n_vals, double);
    BFT_REALLOC(writer->time_steps, n_vals, int);

    writer->time_values[n_vals - 1] = time_value;
    writer->time_steps[n_vals - 1] = time_step;
  }

  return;
}

/*----------------------------------------------------------------------------
 * Indicate if a elements of a given type in a mesh associated to a given
 * MED file writer need to be tesselated.
 *
 * parameters:
 *   this_writer  <-- pointer to associated writer
 *   mesh         <-- pointer to nodal mesh structure that should be written
 *   element_type <-- element type we are interested in
 *
 * returns:
 *   1 if tesselation of the given element type is needed, 0 otherwise
 *----------------------------------------------------------------------------*/

int
fvm_to_med_needs_tesselation(void               *this_writer,
                             const fvm_nodal_t  *mesh,
                             fvm_element_t       element_type)
{
  int  i;
  int  retval = 0;
  fvm_to_med_writer_t  *writer = (fvm_to_med_writer_t *)this_writer;

  /* Initial count and allocation */

  if (   (   element_type == FVM_FACE_POLY
          && writer->divide_polygons == true)
      || (   element_type == FVM_CELL_POLY
          && writer->divide_polyhedra == true)) {

    for (i = 0 ; i < mesh->n_sections ; i++) {

      const fvm_nodal_section_t  *section = mesh->sections[i];

      if (section->type == element_type)
        retval = 1;
    }

  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Indicate a given mesh is present in a MED file.
 *
 * This does not do any verification that the mesh is indeed present,
 * so this should be ensured before. The writer info is simply updated
 * so that additional fields may be output or updated in an existing file.
 *
 * parameters:
 *   this_writer  <-- pointer to associated writer.
 *   mesh         <-- pointer to nodal mesh structure that should be written.
 *----------------------------------------------------------------------------*/

void
fvm_to_med_map_nodal(void               *this_writer,
                     const fvm_nodal_t  *mesh)
{
  int   i_char, n_chars, med_mesh_num;

  char  med_mesh_name[MED_NAME_SIZE + 1];

  fvm_to_med_writer_t  *writer = (fvm_to_med_writer_t *)this_writer;

  /* Clean mesh->name */

  if (mesh->name == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("Mesh name required to continue.\n"));

  strncpy(med_mesh_name, mesh->name, MED_NAME_SIZE);
  n_chars = strlen(med_mesh_name);

  for (i_char = n_chars + 1; i_char < MED_NAME_SIZE; i_char++)
    med_mesh_name[i_char] = ' ';
  med_mesh_name[MED_NAME_SIZE] = '\0';

  /* Get med_mesh_num */

  med_mesh_num = _get_med_mesh_num(writer,
                                   med_mesh_name);

  if (med_mesh_num == 0)
    med_mesh_num = _add_med_mesh(writer,
                                 med_mesh_name,
                                 mesh);
}

/*----------------------------------------------------------------------------
 * Write nodal mesh to a MED file
 *
 * parameters:
 *   this_writer  <-- pointer to associated writer.
 *   mesh         <-- pointer to nodal mesh structure that should be written.
 *----------------------------------------------------------------------------*/

void
fvm_to_med_export_nodal(void               *this_writer,
                        const fvm_nodal_t  *mesh)
{
  int   i_char, n_chars, med_mesh_num;
  med_geometry_type med_type;

  cs_gnum_t   global_connect_slice_size = 0;

  char  med_mesh_name[MED_NAME_SIZE + 1];

  med_int  *med_export_connect = NULL;
  fvm_to_med_mesh_t  *med_mesh = NULL;
  fvm_writer_section_t  *export_list = NULL;

  fvm_to_med_writer_t  *writer = (fvm_to_med_writer_t *)this_writer;

  const fvm_writer_section_t  *export_sections = NULL;
  const int  n_ranks = writer->n_ranks;

  /* Re-open MED file */

  if (writer->is_open == false)
    _med_file_open(writer, MED_ACC_RDWR);

  /* Get med_mesh structure */
  /*------------------------*/

  /* Clean mesh->name */

  if (mesh->name == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("Mesh name required to continue.\n"));

  strncpy(med_mesh_name, mesh->name, MED_NAME_SIZE);
  n_chars = strlen(med_mesh_name);

  for (i_char = n_chars + 1; i_char < MED_NAME_SIZE; i_char++)
    med_mesh_name[i_char] = ' ';
  med_mesh_name[MED_NAME_SIZE] = '\0';

  /* Get med_mesh_num */

  med_mesh_num = _get_med_mesh_num(writer,
                                   med_mesh_name);

  if (med_mesh_num == 0)
    med_mesh_num = _add_med_mesh(writer,
                                 med_mesh_name,
                                 mesh);

  med_mesh = writer->med_meshes[med_mesh_num - 1];

  /*---------------------------*/
  /* Export vertex coordinates */
  /*---------------------------*/

#if defined(HAVE_MPI)
  if (n_ranks > 1)
    _export_vertex_coords_g(mesh,
                            med_mesh,
                            writer);
#endif /* HAVE_MPI */

  if (n_ranks == 1)
    _export_vertex_coords_l(mesh,
                            med_mesh,
                            writer);

  /*---------------------*/
  /* Export connectivity */
  /*---------------------*/

  /* Build list of sections that are used here, in order of output */

  export_list = fvm_writer_export_list(mesh,
                                       0,
                                       true,
                                       false,
                                       writer->discard_polygons,
                                       writer->discard_polyhedra,
                                       writer->divide_polygons,
                                       writer->divide_polyhedra);

  if (export_list == NULL)
    _export_point_elements(writer, mesh,med_mesh);

  /* Compute connectivity buffer size */

  export_sections = export_list;

  if (export_sections != NULL)
    global_connect_slice_size = _get_connect_buffer_size(writer,
                                                         export_sections);

  if (n_ranks == 1)
    BFT_MALLOC(med_export_connect, global_connect_slice_size, med_int);

  /* Loop on MED element types */
  /*---------------------------*/

  while (export_sections != NULL) {

    med_type = _get_med_elt_type(export_sections->type);

#if defined(HAVE_MPI)
    if (n_ranks > 1) { /* Parallel treatment */

      if (med_type == MED_POLYGON) {

        if (writer->discard_polygons == false)
          export_sections =
            _export_nodal_polygons_g(export_sections,
                                     writer,
                                     mesh,
                                     med_mesh);
        else
          /* discard_polygons == true */
          break;

      }
      else if (med_type == MED_POLYHEDRON) {

        if (writer->discard_polyhedra == false)
          export_sections =
            _export_nodal_polyhedra_g(export_sections,
                                      writer,
                                      mesh,
                                      med_mesh);
        else
          /* discard_polyhedra == true */
          break;

      }
      else
        export_sections =
          _export_connect_g(export_sections,
                            writer,
                            mesh,
                            med_mesh);

    }
#endif /* HAVE_MPI */

    if (n_ranks == 1) { /* Serial treatment */

      if (med_type == MED_POLYGON) {

        if (writer->discard_polygons == false)
          export_sections =
            _export_nodal_polygons_l(export_sections,
                                     writer,
                                     med_mesh,
                                     med_export_connect);
        else
          /* discard_polygons == true */
          break;

      }
      else if (med_type == MED_POLYHEDRON) {

        if (writer->discard_polyhedra == false)
          export_sections =
            _export_nodal_polyhedra_l(export_sections,
                                      writer,
                                      med_mesh,
                                      med_export_connect);
        else
          /* discard_polyhedra == true */
          break;

      }
      else
        export_sections =
          _export_connect_l(export_sections,
                            writer,
                            med_mesh,
                            med_export_connect);

    } /* end of serial treatment */

  } /* End of loop on med_types */

  /* Free buffers */

  if (med_export_connect != NULL)
    BFT_FREE(med_export_connect);

  if (export_list != NULL)
    BFT_FREE(export_list);

  /* Close MED file (to force its update) */

  if (writer->is_open == true)
    _med_file_close(writer);
}

/*----------------------------------------------------------------------------
 * Write field associated with a nodal mesh to a MED file.
 *
 * Assigning a negative value to the time step indicates a time-independent
 * field (in which case the time_value argument is unused).
 *
 * parameters:
 *   this_writer      <-- pointer to associated writer
 *   mesh             <-- pointer to associated nodal mesh structure
 *   name             <-- variable name
 *   location         <-- variable definition location (nodes or elements)
 *   dimension        <-- variable dimension (0: constant, 1: scalar,
 *                        3: vector, 6: sym. tensor, 9: asym. tensor)
 *   interlace        <-- indicates if variable in memory is interlaced
 *   n_parent_lists   <-- indicates if variable values are to be obtained
 *                        directly through the local entity index (when 0) or
 *                        through the parent entity numbers (when 1 or more)
 *   parent_num_shift <-- parent number to value array index shifts;
 *                        size: n_parent_lists
 *   datatype         <-- indicates the data type of (source) field values
 *   time_step        <-- number of the current time step
 *   time_value       <-- associated time value
 *   field_values     <-- array of associated field value arrays
 *----------------------------------------------------------------------------*/

void
fvm_to_med_export_field(void                            *this_writer,
                        const fvm_nodal_t               *mesh,
                        const char                      *name,
                        const fvm_writer_var_loc_t       location,
                        const int                        dimension,
                        const cs_interlace_t             interlace,
                        const int                        n_parent_lists,
                        const cs_lnum_t                  parent_num_shift[],
                        const cs_datatype_t              datatype,
                        const int                        time_step,
                        const double                     time_value,
                        const void                *const field_values[])
{
  int   i_char, n_chars, data_sizeof;
  int   med_mesh_num;
  char  med_mesh_name[MED_NAME_SIZE + 1];
  char  med_fieldname[MED_NAME_SIZE + 1];

  cs_datatype_t  datatype_convert = CS_DATATYPE_NULL;
  med_field_type  datatype_med = MED_FLOAT64;

  fvm_to_med_mesh_t *med_mesh = NULL;
  fvm_writer_section_t  *export_list = NULL;
  fvm_to_med_writer_t  *writer = (fvm_to_med_writer_t *)this_writer;

  const fvm_writer_section_t  *export_sections = NULL;
  fvm_writer_field_helper_t  *helper = NULL;

  /* Re-open MED file */

  if (writer->is_open == false)
    _med_file_open(writer, MED_ACC_RDWR);

  /* Adapt dimension */

  if (dimension != 1 && dimension != 3 && dimension != 6 && dimension != 9)
    bft_error(__FILE__, __LINE__, 0,
              _("Data of dimension %d not handled"), dimension);

  /* Get med_mesh_name */

  strncpy(med_mesh_name, mesh->name, MED_NAME_SIZE);
  n_chars = strlen(med_mesh_name);

  for (i_char = n_chars + 1; i_char < MED_NAME_SIZE; i_char++)
    med_mesh_name[i_char] = ' ';

  med_mesh_name[MED_NAME_SIZE] = '\0';

  /* Get MED mesh structure */
  /*------------------------*/

  med_mesh_num = _get_med_mesh_num(writer,
                                   med_mesh_name);

  if (med_mesh_num == 0)
    med_mesh_num = _add_med_mesh(writer,
                                 med_mesh_name,
                                 mesh);

  med_mesh = writer->med_meshes[med_mesh_num - 1];

  /* Adapt cs_datatype and find corresponding MED datatype */

  _get_datatypes(datatype,
                 &datatype_convert,
                 &datatype_med,
                 &data_sizeof);

  /* Create MED field if necessary. Always return MED fieldname */

  _get_med_fieldname(writer,
                     med_mesh_name,
                     name,
                     datatype_med,
                     dimension,
                     med_fieldname);

  /* Build list of sections that are used here, in order of output */
  /*---------------------------------------------------------------*/

  export_list = fvm_writer_export_list(mesh,
                                       fvm_nodal_get_max_entity_dim(mesh),
                                       true,
                                       false,
                                       writer->discard_polygons,
                                       writer->discard_polyhedra,
                                       writer->divide_polygons,
                                       writer->divide_polyhedra);

  /* Build writer helper */
  /*---------------------*/

  helper = fvm_writer_field_helper_create(mesh,
                                          export_list,
                                          dimension,
                                          CS_INTERLACE,
                                          datatype_convert,
                                          location);

#if defined(HAVE_MPI)

  if (writer->n_ranks > 1)
    fvm_writer_field_helper_init_g(helper,
                                   writer->min_rank_step,
                                   writer->min_block_size,
                                   writer->comm);

#endif

  /* Export field */
  /*--------------*/

  if (location == FVM_WRITER_PER_NODE)
    _export_field_values_n(mesh,
                           writer,
                           helper,
                           dimension,
                           interlace,
                           n_parent_lists,
                           parent_num_shift,
                           datatype,
                           field_values,
                           med_mesh->name,
                           med_fieldname,
                           time_step,
                           time_value);

  else if (location == FVM_WRITER_PER_ELEMENT) {

    if (export_list == NULL)
      bft_error(__FILE__, __LINE__, 0,
                _("MED must have entities.\n"
                  "Mesh: \"%s\"\n"
                  "Writer: \"%s\"\n"),
                med_mesh->name, writer->name);

    /* Loop on MED element types */

    export_sections = export_list;

    while (export_sections != NULL) {

      export_sections = _export_field_values_e(export_sections,
                                               writer,
                                               helper,
                                               dimension,
                                               interlace,
                                               n_parent_lists,
                                               parent_num_shift,
                                               datatype,
                                               field_values,
                                               med_mesh->name,
                                               med_fieldname,
                                               time_step,
                                               time_value);

    } /* End of loop on MED element types */

  }
  else
    bft_error(__FILE__, __LINE__, 0,
              "fvm_to_med_export_field(): field location not managed.\n"
              "Associated writer: \"%s\"\n"
              "Associated med_mesh: \"%s\"\n"
              "Associated fieldname: \"%s\"\n"
              "Associated location: %i\n",
              writer->name, med_mesh_name, med_fieldname, location);

  /* Free helper structures */

  fvm_writer_field_helper_destroy(&helper);

  BFT_FREE(export_list);

  /* Close MED file (to force its update) */
  /*--------------------------------------*/

  _med_file_close(writer);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* HAVE_MED */
