/*============================================================================
 * Write a nodal representation associated with a mesh and associated
 * variables to MEDCoupling objects
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

#if defined(HAVE_MEDCOUPLING)

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

#include <MEDCoupling_version.h>

#include <MEDCouplingUMesh.hxx>
#include <MEDCouplingField.hxx>
#include <MEDCouplingFieldDouble.hxx>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "fvm_defs.h"
#include "fvm_convert_array.h"
#include "fvm_io_num.h"
#include "fvm_nodal.h"
#include "fvm_nodal_priv.h"
#include "fvm_writer_priv.h"

#include "cs_block_dist.h"
#include "cs_file.h"
#include "cs_parall.h"
#include "cs_part_to_block.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvm_to_medcoupling.h"

/*----------------------------------------------------------------------------*/

using namespace MEDCoupling;

/*============================================================================
 * Local Type Definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * MED field structure
 *----------------------------------------------------------------------------*/

typedef struct {

  int                       mesh_id;       /* Associated mesh structure id */

  int                       dim;           /* Field dimension */
  TypeOfTimeDiscretization  td;            /* NO_TIME, ONE_TIME, LINEAR_TIME,
                                              or CONST_ON_TIME_INTERVAL */
  MEDCouplingFieldDouble   *f;             /* Pointer to MED writer field */

} fvm_medcoupling_field_t;

/*----------------------------------------------------------------------------
 * MEDCoupling writer/reader structure
 *----------------------------------------------------------------------------*/

typedef struct {

  char                      *name;           /* Writer name */

  int                        rank;           /* Rank of current process
                                                in communicator */
  int                        n_ranks;        /* Size of communicator */

  int                        n_med_meshes;   /* Number of MED meshes */
  MEDCouplingUMesh         **med_meshes;     /* Array of pointers to MED mesh
                                                structure */
  fvm_writer_time_dep_t      time_dependency; /* Mesh time dependency */

  int                        n_fields;       /* Number of fields */
  fvm_medcoupling_field_t  **fields;         /* Array of field helper
                                                structures */

  int                        n_time_steps;   /* Number of meshe time steps */
  int                       *time_steps;     /* Array of meshe time steps */
  double                    *time_values;    /* Array of meshe time values */

#if defined(HAVE_MPI)
  MPI_Comm                   comm;           /* Associated MPI communicator */
#endif

} fvm_to_medcoupling_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Return the MEDCoupling mesh id associated with a given mesh name,
 * or -1 if no association found.
 *
 * parameters:
 *   writer    <-- writer structure
 *   mesh_name <-- mesh name
 *
 * returns:
 *    MEDCouping mesh id, or -1 if MEDCoupling mesh name is not associated
 *    with this writer structure in FVM
 *----------------------------------------------------------------------------*/

static int
_get_medcoupling_mesh_id(fvm_to_medcoupling_t  *writer,
                         const char            *mesh_name)
{
  int i;
  int retval = -1;

  assert(writer != NULL);

  for (i = 0; i < writer->n_med_meshes; i++) {
    if (strcmp(mesh_name, writer->med_meshes[i]->getName().c_str()) == 0)
      break;
  }

  if (i < writer->n_med_meshes)
    retval = i;

  return retval;
}

/*----------------------------------------------------------------------------
 * Create a MEDCoupling mesh structure.
 *
 * parameters:
 *   writer    <-- MEDCoupling structure.
 *   mesh      <-- FVM mesh  structure.
 *
 * returns:
 *   MEDCoupling mesh id.
 *----------------------------------------------------------------------------*/

static int
_add_medcoupling_mesh(fvm_to_medcoupling_t  *writer,
                      const fvm_nodal_t     *mesh)
{
  int id;

  assert(writer != NULL);

  /* Add a new MEDCoupling mesh structure */

  MEDCouplingUMesh *m = NULL;

  m = MEDCouplingUMesh::New();

  m->setName(mesh->name);
  m->setDescription("Generated by code_saturne/FVM.");
  m->setTimeUnit("s");
  m->setMeshDimension(fvm_nodal_get_max_entity_dim(mesh));

  writer->n_med_meshes += 1;
  id = writer->n_med_meshes - 1;

  BFT_REALLOC(writer->med_meshes, writer->n_med_meshes, MEDCouplingUMesh *);

  writer->med_meshes[id] = m;

  return (writer->n_med_meshes - 1);
}

/*----------------------------------------------------------------------------
 * Define MED geometrical element type according to FVM element type
 *
 * parameters:
 *   fvm_elt_type <-- pointer to fvm element type.
 *
 * return:
 *   med geometrical element type.
 *----------------------------------------------------------------------------*/

static INTERP_KERNEL::NormalizedCellType
_get_norm_elt_type(const fvm_element_t fvm_elt_type)
{
  INTERP_KERNEL::NormalizedCellType  norm_elt_type;

  switch (fvm_elt_type) {

  case FVM_EDGE:
    norm_elt_type = INTERP_KERNEL::NORM_SEG2;
    break;

  case FVM_FACE_TRIA:
    norm_elt_type = INTERP_KERNEL::NORM_TRI3;
    break;

  case FVM_FACE_QUAD:
    norm_elt_type = INTERP_KERNEL::NORM_QUAD4;
    break;

  case FVM_FACE_POLY:
    norm_elt_type = INTERP_KERNEL::NORM_POLYGON;
    break;

  case FVM_CELL_TETRA:
    norm_elt_type = INTERP_KERNEL::NORM_TETRA4;
    break;

  case FVM_CELL_PYRAM:
    norm_elt_type = INTERP_KERNEL::NORM_PYRA5;
    break;

  case FVM_CELL_PRISM:
    norm_elt_type = INTERP_KERNEL::NORM_PENTA6;
    break;

  case FVM_CELL_HEXA:
    norm_elt_type = INTERP_KERNEL::NORM_HEXA8;
    break;

  case FVM_CELL_POLY:
    norm_elt_type = INTERP_KERNEL::NORM_POLYHED;
    break;

  default:
    norm_elt_type = INTERP_KERNEL::NORM_ERROR;
    bft_error(__FILE__, __LINE__, 0,
              "_get_norm_elt_type(): "
              "No association with INTERP_KERNEL element type has been found\n"
              "FVM element type: \"%i\"\n",
              (int)fvm_elt_type);

  } /* End of switch on element type */

  return norm_elt_type;
}

/*----------------------------------------------------------------------------
 * Get vertex order to describe MED element type.
 *
 * parameters:
 *   norm_elt_type  <-- MED element type.
 *   vertex_order  --> Pointer to vertex order array (0 to n-1).
 *
 *----------------------------------------------------------------------------*/

static void
_get_vertex_order(INTERP_KERNEL::NormalizedCellType   norm_elt_type,
                  int                                *vertex_order)
{
  switch(norm_elt_type) {

  case INTERP_KERNEL::NORM_SEG2:
    vertex_order[0] = 0;
    vertex_order[1] = 1;
    break;

  case INTERP_KERNEL::NORM_TRI3:
    vertex_order[0] = 0;
    vertex_order[1] = 2;
    vertex_order[2] = 1;
    break;

  case INTERP_KERNEL::NORM_QUAD4:
    vertex_order[0] = 0;
    vertex_order[1] = 3;
    vertex_order[2] = 2;
    vertex_order[3] = 1;
    break;

  case INTERP_KERNEL::NORM_TETRA4:
    vertex_order[0] = 0;
    vertex_order[1] = 2;
    vertex_order[2] = 1;
    vertex_order[3] = 3;
    break;

  case INTERP_KERNEL::NORM_PYRA5:
    vertex_order[0] = 0;
    vertex_order[1] = 3;
    vertex_order[2] = 2;
    vertex_order[3] = 1;
    vertex_order[4] = 4;
    break;

  case INTERP_KERNEL::NORM_PENTA6:
    vertex_order[0] = 0;
    vertex_order[1] = 2;
    vertex_order[2] = 1;
    vertex_order[3] = 3;
    vertex_order[4] = 5;
    vertex_order[5] = 4;
    break;

  case INTERP_KERNEL::NORM_HEXA8:
    vertex_order[0] = 0;
    vertex_order[1] = 3;
    vertex_order[2] = 2;
    vertex_order[3] = 1;
    vertex_order[4] = 4;
    vertex_order[5] = 7;
    vertex_order[6] = 6;
    vertex_order[7] = 5;
    break;

  case INTERP_KERNEL::NORM_POLYGON:
    vertex_order[0] = -1;
    break;

  case INTERP_KERNEL::NORM_POLYHED:
    vertex_order[0] = -1;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "_get_vertex_order(): No associated MED element type known\n"
              "INTERP_KERNEL element type: \"%i\"\n",
              (int)norm_elt_type);
  }

  return;
}

/*----------------------------------------------------------------------------
 * Return the MEDCoupling field id associated with given mesh and field names,
 * or -1 if no association found.
 *
 * parameters:
 *   writer    <-- MED writer structure.
 *   fieldname <-- input fieldname.
 *   mesh_id   <-- id of associated mesh in structure.
 *   dim       <-- number of field components.
 *   location  <-- mesh location (cells or vertices).
 *   type      <-- field type (cells, nodes)
 *   td        <-- time discretization type
 *
 * returns
 *   field id in writer structure, or -1 if not found
 *----------------------------------------------------------------------------*/

static int
_get_medcoupling_field_id(fvm_to_medcoupling_t      *writer,
                          const char                *fieldname,
                          int                        mesh_id,
                          int                        dim,
                          fvm_writer_var_loc_t       location,
                          TypeOfTimeDiscretization   td)
{
  int n_fields, f_id;

  TypeOfField type = (location == FVM_WRITER_PER_NODE) ? ON_NODES : ON_CELLS;

  MEDCouplingFieldDouble  *f = NULL;

  /* Loop on fields to know if field has already been created */

  n_fields = writer->n_fields;

  for (f_id = 0; f_id < n_fields; f_id++) {

    f = (writer->fields[f_id])->f;

    if (   writer->fields[f_id]->mesh_id == mesh_id
        && strcmp(fieldname, f->getName().c_str()) == 0) {

      /* If field exists, check that dimensions and type are compatible */

      fvm_medcoupling_field_t *field = writer->fields[f_id];
      TypeOfField type_ref = f->getTypeOfField();

      if (field->dim != dim)
        bft_error(__FILE__, __LINE__, 0,
                  _("MEDCoupling field \"%s\" already defined for\n"
                    "coupling \"%s\" and mesh \"%s\" with %d components,\n"
                    "but re-defined with %d components."),
                  fieldname, writer->name,
                  writer->med_meshes[field->mesh_id]->getName().c_str(),
                  field->dim, dim);

      else if (type_ref != type)
        bft_error(__FILE__, __LINE__, 0,
                  _("MEDCoupling field \"%s\" already defined for\n"
                    "coupling \"%s\" and mesh \"%s\" with type %d,\n"
                    "but re-defined with type %d."),
                  fieldname, writer->name,
                  writer->med_meshes[field->mesh_id]->getName().c_str(),
                  (int)type_ref, type);

      else if ((writer->fields[f_id])->td != td)
        bft_error(__FILE__, __LINE__, 0,
                  _("MEDCoupling field \"%s\" already defined for writer\n"
                    "\"%s\" and mesh \"%s\" with time discretization %d,\n"
                    "but re-defined with time discretization %d."),
                  fieldname, writer->name,
#if (MEDCOUPLING_VERSION  >= 0x070300)
                  writer->med_meshes[field->mesh_id]->getName().c_str(),
#else
                  writer->med_meshes[field->mesh_id]->getName(),
#endif
                  (int)(field->td), (int)td);

      /* return id of field if compatible */

      else
        return f_id;

      break;
    }
  }

  return -1;
}

/*----------------------------------------------------------------------------
 * Create a MEDCoupling field structure.
 *
 * parameters:
 *   writer  <-- MED writer structure.
 *   fieldname <-- input fieldname.
 *   mesh_id   <-- id of associated mesh in structure.
 *   dim       <-- number of field components.
 *   location  <-- mesh location (cells or vertices).
 *   type      <-- field type (cells, nodes)
 *   td        <-- time discretization type
 *
 * returns
 *   field id in writer structure
 *----------------------------------------------------------------------------*/

static int
_add_medcoupling_field(fvm_to_medcoupling_t      *writer,
                       const char                *fieldname,
                       int                        mesh_id,
                       int                        dim,
                       fvm_writer_var_loc_t       location,
                       TypeOfTimeDiscretization   td)
{
  TypeOfField type = (location == FVM_WRITER_PER_NODE) ? ON_NODES : ON_CELLS;

  int f_id = -1;

  /* Prepare writer structure */

  f_id = writer->n_fields;

  BFT_REALLOC(writer->fields,
              writer->n_fields + 1,
              fvm_medcoupling_field_t *);

  /* Build ParaFIELD object if required */

  MEDCouplingFieldDouble  *f = NULL;

  if (writer->med_meshes[mesh_id] != NULL) {

    f = MEDCouplingFieldDouble::New(type, td);

    f->setName(fieldname);

    /* Assign array to field (filled later) */

    int n_locs = 0;
    DataArrayDouble *array = DataArrayDouble::New();

    if (location == FVM_WRITER_PER_NODE)
      n_locs = writer->med_meshes[mesh_id]->getNumberOfNodes();
    else if (location == FVM_WRITER_PER_ELEMENT)
      n_locs = writer->med_meshes[mesh_id]->getNumberOfCells();

    array->alloc(n_locs, dim);
    f->setMesh(writer->med_meshes[mesh_id]);
    f->setArray(array);
    array->decrRef();

  }

  /* Update writer structure */

  f_id = writer->n_fields;

  BFT_REALLOC(writer->fields,
              writer->n_fields + 1,
              fvm_medcoupling_field_t *);

  BFT_MALLOC(writer->fields[f_id], 1, fvm_medcoupling_field_t);

  writer->fields[f_id]->mesh_id = mesh_id;
  writer->fields[f_id]->td = td;
  writer->fields[f_id]->dim = dim;

  writer->fields[f_id]->f = f;

  writer->n_fields++;

  return f_id;
}

/*----------------------------------------------------------------------------
 * Write vertex coordinates to a MEDCoupling object in serial mode
 *
 * parameters:
 *   mesh        <-- pointer to nodal mesh structure
 *   med_mesh    <-- pointer to MEDCouuplingUMesh object (NULL on ranks > 0)
 *----------------------------------------------------------------------------*/

static void
_export_vertex_coords_l(const fvm_nodal_t           *mesh,
                        MEDCouplingUMesh            *med_mesh)
{
  cs_lnum_t   i, j;
  size_t stride;

  DataArrayDouble  *Coords = NULL;
  double  *block_coords = NULL;

  const double  *vertex_coords = mesh->vertex_coords;
  const cs_lnum_t  n_vertices = mesh->n_vertices;

  /* Vertex coordinates */
  /*--------------------*/

  stride = (size_t)(mesh->dim);

  Coords = DataArrayDouble::New();
  Coords->alloc(n_vertices, 3);
  block_coords = Coords->getPointer();

  if (mesh->parent_vertex_num != NULL || mesh->dim < 3) {

    if (mesh->parent_vertex_num != NULL) {
      const cs_lnum_t  *parent_vertex_num = mesh->parent_vertex_num;
      for (i = 0; i < n_vertices; i++) {
        for (j = 0; j < mesh->dim; j++)
          block_coords[i*3 + j]
            = vertex_coords[(parent_vertex_num[i]-1)*stride + j];
        for (; j < 3; j++)
          block_coords[i*3 + j] = 0.;
      }
    }
    else {
      for (i = 0; i < n_vertices; i++) {
        for (j = 0; j < mesh->dim; j++)
          block_coords[i*3 + j] = vertex_coords[i*stride + j];
        for (; j < 3; j++)
          block_coords[i*3 + j] = 0.;
      }
    }

  }
  else
    memcpy(block_coords, vertex_coords, 3*mesh->n_vertices*sizeof(double));

  med_mesh->setCoords(Coords);
  Coords->decrRef();
}

/*----------------------------------------------------------------------------
 * Write strided connectivity block to a MEDCoupling mesh
 *
 * The connectivity on input may use 1 to n numbering, so it is shifted
 * by -1 here.
 *
 * parameters:
 *   type     <-- FVM element type
 *   n_elts   <-- number of elements in block
 *   connect  <-- connectivity array
 *   med_mesh <-> pointer to MEDCouuplingUMesh object (NULL on ranks > 0)
 *----------------------------------------------------------------------------*/

static void
_write_connect_block(fvm_element_t      type,
                     cs_lnum_t          n_elts,
                     const cs_lnum_t    connect[],
                     MEDCouplingUMesh  *med_mesh)
{
  int vertex_order[8];
  mcIdType elt_buf[8];
  cs_lnum_t  i;
  int  j;

  const mcIdType  stride = fvm_nodal_n_vertices_element[type];
  INTERP_KERNEL::NormalizedCellType med_type = _get_norm_elt_type(type);

  _get_vertex_order(med_type, vertex_order);

  assert(med_mesh != NULL);

  for (i = 0; i < n_elts; i++) {
    for (j = 0; j < stride; j++)
      elt_buf[j] = connect[i*stride + vertex_order[j]] - 1;
    med_mesh->insertNextCell(med_type, stride, elt_buf);
  }
}

/*----------------------------------------------------------------------------
 * Write polyhedra from a nodal mesh to a MEDCoupling object in serial mode
 *
 * parameters:
 *   export_section <-- pointer to MEDCoupling section helper structure
 *   med_mesh       <-> MEDCouuplingUMesh object (NULL on ranks > 0)
 *----------------------------------------------------------------------------*/

static void
_export_nodal_polyhedra_l(const fvm_nodal_section_t  *section,
                          MEDCouplingUMesh           *med_mesh)
{
  int  face_sgn;
  cs_lnum_t  i, j, k, l;

  cs_lnum_t  face_length, face_id;

  int elt_buf_size = 8;
  mcIdType *elt_buf = NULL;

  BFT_MALLOC(elt_buf, elt_buf_size, mcIdType);

  /* Write cell/vertex connectivity */
  /*--------------------------------*/

  for (i = 0; i < section->n_elements; i++) {

    mcIdType  m = 0;

    /* Loop on cell faces */

    for (j = section->face_index[i];
         j < section->face_index[i+1];
         j++) {

      /* Print face vertex numbers */

      if (section->face_num[j] > 0) {
        face_id = section->face_num[j] - 1;
        face_sgn = 1;
      }
      else {
        face_id = -section->face_num[j] - 1;
        face_sgn = -1;
      }

      face_length = (  section->vertex_index[face_id+1]
                     - section->vertex_index[face_id]);

      while (m + face_length + 1 > elt_buf_size) {
        elt_buf_size *= 2;
        BFT_REALLOC(elt_buf, elt_buf_size, mcIdType);
      }

      if (j >  section->face_index[i])
        elt_buf[m++] = -1;

      for (k = 0; k < face_length; k++) {
        l =    section->vertex_index[face_id]
            + (face_length + (k*face_sgn))%face_length;
        elt_buf[m++] = section->vertex_num[l] - 1;
      }

    } /* End of loop on cell faces */

    med_mesh->insertNextCell(INTERP_KERNEL::NORM_POLYHED, m, elt_buf);

  } /* End of loop on polyhedral cells */

  BFT_FREE(elt_buf);
}

/*----------------------------------------------------------------------------
 * Write polygons from a nodal mesh to a text object in serial mode
 *
 * parameters:
 *   export_section <-- pointer to MEDCoupling section helper structure
 *   med_mesh       <-> MEDCouuplingUMesh object
*----------------------------------------------------------------------------*/

static void
_export_nodal_polygons_l(const fvm_nodal_section_t  *section,
                         MEDCouplingUMesh           *med_mesh)


{
  cs_lnum_t   i, j;

  int elt_buf_size = 8;
  mcIdType *elt_buf = NULL;

  BFT_MALLOC(elt_buf, elt_buf_size, mcIdType);

  /* Loop on all polygonal faces */
  /*-----------------------------*/

  for (i = 0; i < section->n_elements; i++) {

    mcIdType k = 0;

    int face_size = section->vertex_index[i+1] - section->vertex_index[i];
    while (elt_buf_size < face_size) {
      elt_buf_size *= 2;
      BFT_REALLOC(elt_buf, elt_buf_size, mcIdType);
    }

    for (j = section->vertex_index[i];
         j < section->vertex_index[i+1];
         j++)
      elt_buf[k++] = section->vertex_num[j] - 1;

    med_mesh->insertNextCell(INTERP_KERNEL::NORM_POLYGON, k, elt_buf);

  } /* End of loop on polygonal faces */

  BFT_FREE(elt_buf);
}

/*----------------------------------------------------------------------------
 * Write field values associated with nodal values of a nodal mesh to
 * a MEDCoupling object in serial mode.
 *
 * Output fields ar either scalar or 3d vectors or scalars, and are
 * non interlaced. Input arrays may be less than 2d, in which case the z
 * values are set to 0, and may be interlaced or not.
 *
 * parameters:
 *   mesh             <-- pointer to nodal mesh structure
 *   dim              <-- field dimension
 *   interlace        <-- indicates if field in memory is interlaced
 *   n_parent_lists   <-- indicates if field values are to be obtained
 *                        directly through the local entity index (when 0) or
 *                        through the parent entity numbers (when 1 or more)
 *   parent_num_shift <-- parent list to common number index shifts;
 *                        size: n_parent_lists
 *   datatype         <-- input data type (output is real)
 *   field_values     <-- array of associated field value arrays
 *   f                <-- associated MEDCouplingFieldDouble object
 *----------------------------------------------------------------------------*/

static void
_export_field_values_n(const fvm_nodal_t           *mesh,
                       int                          dim,
                       cs_interlace_t               interlace,
                       int                          n_parent_lists,
                       const cs_lnum_t              parent_num_shift[],
                       cs_datatype_t                datatype,
                       const void            *const field_values[],
                       MEDCouplingFieldDouble      *f)
{
  assert(f != NULL);

  double  *values = f->getArray()->getPointer();

  fvm_convert_array(dim,
                    0,
                    dim, /* stride */
                    0, /* start_id */
                    mesh->n_vertices, /* end_id */
                    interlace,
                    datatype,
                    CS_DOUBLE,
                    n_parent_lists,
                    parent_num_shift,
                    mesh->parent_vertex_num,
                    field_values,
                    values);
}

/*----------------------------------------------------------------------------
 * Write field values associated with element values of a nodal mesh to
 * a MEDCoupling object.
 *
 * Output fields are non interlaced. Input arrays may be interlaced or not.
 *
 * parameters:
 *   mesh             <-- pointer to nodal mesh structure
 *   dim              <-- field dimension
 *   interlace        <-- indicates if field in memory is interlaced
 *   n_parent_lists   <-- indicates if field values are to be obtained
 *                        directly through the local entity index (when 0) or
 *                        through the parent entity numbers (when 1 or more)
 *   parent_num_shift <-- parent list to common number index shifts;
 *                        size: n_parent_lists
 *   datatype         <-- indicates the data type of (source) field values
 *   field_values     <-- array of associated field value arrays
 *   f                <-- associated object handle
 *----------------------------------------------------------------------------*/

static void
_export_field_values_e(const fvm_nodal_t               *mesh,
                       int                              dim,
                       cs_interlace_t                   interlace,
                       int                              n_parent_lists,
                       const cs_lnum_t                  parent_num_shift[],
                       cs_datatype_t                    datatype,
                       const void                *const field_values[],
                       MEDCouplingFieldDouble          *f)
{
  int  section_id;

  double  *values = NULL;

  if (f != NULL)
    values = f->getArray()->getPointer();

  /* Distribute partition to block values */

  cs_lnum_t start_id = 0;
  cs_lnum_t src_shift = 0;

  /* loop on sections which should be appended */

  const int  elt_dim = fvm_nodal_get_max_entity_dim(mesh);

  for (section_id = 0; section_id < mesh->n_sections; section_id++) {

    const fvm_nodal_section_t  *section = mesh->sections[section_id];

    if (section->entity_dim < elt_dim)
      continue;

    fvm_convert_array(dim,
                      0,
                      dim,
                      src_shift,
                      section->n_elements + src_shift,
                      interlace,
                      datatype,
                      CS_DOUBLE,
                      n_parent_lists,
                      parent_num_shift,
                      section->parent_element_num,
                      field_values,
                      values + start_id);

    start_id += section->n_elements*dim;
    if (n_parent_lists == 0)
      src_shift += section->n_elements;

  }
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif

/*----------------------------------------------------------------------------
 * Initialize FVM to MEDCoupling object writer.
 *
 * No options are available for this format.
 *
 * parameters:
 *   name           <-- base output case name.
 *   options        <-- whitespace separated, lowercase options list
 *   time_dependecy <-- indicates if and how meshes will change with time
 *   comm           <-- associated MPI communicator.
 *
 * returns:
 *   pointer to opaque MEDCoupling writer structure.
 *----------------------------------------------------------------------------*/

#if defined(HAVE_MPI)
void *
fvm_to_medcoupling_init_writer(const char             *name,
                               const char             *path,
                               const char             *options,
                               fvm_writer_time_dep_t   time_dependency,
                               MPI_Comm                comm)
#else
void *
fvm_to_medcoupling_init_writer(const char             *name,
                               const char             *path,
                               const char             *options,
                               fvm_writer_time_dep_t   time_dependency)
#endif
{
  CS_UNUSED(path);
  CS_UNUSED(options);

  fvm_to_medcoupling_t  *writer = NULL;

  /* Initialize writer */

  BFT_MALLOC(writer, 1, fvm_to_medcoupling_t);

  writer->rank = 0;
  writer->n_ranks = 1;

  writer->n_med_meshes = 0;
  writer->n_fields  = 0;
  writer->med_meshes   = NULL;
  writer->fields = NULL;

  writer->n_time_steps   = 0;
  writer->time_steps     = NULL;
  writer->time_values    = NULL;
  writer->time_dependency = time_dependency;

  /* Writer name */

  if (name != NULL) {
    BFT_MALLOC(writer->name, strlen(name) + 1, char);
    strcpy(writer->name, name);
  }
  else {
    const char _name[] = "MEDCoupling writer";
    BFT_MALLOC(writer->name, strlen(_name) + 1, char);
    strcpy(writer->name, _name);
  }

  /* Parallel parameters */

#if defined(HAVE_MPI)
  {
    int mpi_flag;
    MPI_Initialized(&mpi_flag);

    if (mpi_flag && comm != MPI_COMM_NULL) {
      writer->comm = comm;
      MPI_Comm_rank(writer->comm, &(writer->rank));
      MPI_Comm_size(writer->comm, &(writer->n_ranks));
    }
    else
      writer->comm = MPI_COMM_NULL;
  }
#endif /* defined(HAVE_MPI) */

  return writer;
}

/*----------------------------------------------------------------------------
 * Finalize FVM to MEDCoupling object writer.
 *
 * parameters:
 *   this_writer_p <-- pointer to opaque writer structure.
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

void *
fvm_to_medcoupling_finalize_writer(void  *this_writer_p)
{
  int i;

  fvm_to_medcoupling_t  *writer = (fvm_to_medcoupling_t *)this_writer_p;

  assert(writer != NULL);

  /* Free structures */

  BFT_FREE(writer->name);
  BFT_FREE(writer->time_values);
  BFT_FREE(writer->time_steps);

  /* Free MEDCouplingUMesh and field structures
     (reference counters should go to 0) */

  for (i = 0; i < writer->n_med_meshes; i++)
    writer->med_meshes[i] = NULL; // delete writer->med_meshes[i];
  BFT_FREE(writer->med_meshes);

  for (i = 0; i < writer->n_fields; i++) {
    writer->fields[i]->f = NULL; // delete writer->fields[i]->f;
    BFT_FREE(writer->fields[i]);
  }

  BFT_FREE(writer->fields);

  /* Free fvm_to_medcoupling_t structure */

  BFT_FREE(writer);

  return NULL;
}

/*----------------------------------------------------------------------------
 * Associate new time step with an MEDCoupling geometry.
 *
 * parameters:
 *   this_writer_p <-- pointer to associated writer
 *   time_step     <-- time step number
 *   time_value    <-- time_value number
 *----------------------------------------------------------------------------*/

void
fvm_to_medcoupling_set_mesh_time(void          *this_writer_p,
                                 const int      time_step,
                                 const double   time_value)
{
  CS_UNUSED(time_step);
  CS_UNUSED(time_value);

  fvm_to_medcoupling_t  *w = (fvm_to_medcoupling_t *)this_writer_p;

  for (int i = 0; i < w->n_med_meshes; i++) {
    MEDCouplingUMesh  *m = w->med_meshes[i];
    m->updateTime();
  }
}

/*----------------------------------------------------------------------------
 * Write nodal mesh to a MEDCoupling object
 *
 * parameters:
 *   this_writer_p <-- pointer to associated writer
 *   mesh          <-- pointer to nodal mesh structure that should be written
 *----------------------------------------------------------------------------*/

void
fvm_to_medcoupling_export_nodal(void               *this_writer_p,
                                const fvm_nodal_t  *mesh)
{
  int  mesh_id, section_id;

  cs_lnum_t  n_g_elts = 0;
  fvm_to_medcoupling_t  *this_writer = (fvm_to_medcoupling_t *)this_writer_p;

  const int  elt_dim = fvm_nodal_get_max_entity_dim(mesh);

  /* Initialization */
  /*----------------*/

  /* Get matching mesh */

  mesh_id = _get_medcoupling_mesh_id(this_writer,
                                     mesh->name);

  if (mesh_id < 0)
    mesh_id = _add_medcoupling_mesh(this_writer,
                                    mesh);

  MEDCouplingUMesh  *med_mesh = this_writer->med_meshes[mesh_id];

  /* Vertex coordinates */
  /*--------------------*/

  _export_vertex_coords_l(mesh, med_mesh);

  /* Element connectivity size */
  /*---------------------------*/

  for (section_id = 0; section_id < mesh->n_sections; section_id++) {

    const fvm_nodal_section_t  *section = mesh->sections[section_id];

    if (section->entity_dim < elt_dim)
      continue;

    if (section->global_element_num != NULL)
      n_g_elts += fvm_io_num_get_global_count(section->global_element_num);
    else
      n_g_elts += section->n_elements;

  } /* End of loop on sections */

  if (med_mesh != NULL)
    med_mesh->allocateCells(n_g_elts);

  /* Element connectivity */
  /*----------------------*/

  for (section_id = 0; section_id < mesh->n_sections; section_id++) {

    const fvm_nodal_section_t  *section = mesh->sections[section_id];

    if (section->entity_dim < elt_dim)
      continue;

    /* Output for strided (regular) element types */
    /*--------------------------------------------*/

    if (section->stride > 0)
      _write_connect_block(section->type,
                           section->n_elements,
                           section->vertex_num,
                           med_mesh);

    /* Output for polygons */
    /*---------------------*/

    else if (section->type == FVM_FACE_POLY)
      _export_nodal_polygons_l(section, med_mesh);

    /* Output for polyhedra */
    /*----------------------*/

    else if (section->type == FVM_CELL_POLY)
      _export_nodal_polyhedra_l(section, med_mesh);

  } /* End of loop on sections */

  /* Update mesh object */
  /*--------------------*/

  med_mesh->finishInsertingCells();
}

/*----------------------------------------------------------------------------
 * Write field associated with a nodal mesh to a MEDCoupling object.
 *
 * Assigning a negative value to the time step indicates a time-independent
 * field (in which case the time_value argument is unused).
 *
 * parameters:
 *   this_writer_p    <-- pointer to associated writer
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
fvm_to_medcoupling_export_field(void                  *this_writer_p,
                                const fvm_nodal_t     *mesh,
                                const char            *name,
                                fvm_writer_var_loc_t   location,
                                int                    dimension,
                                cs_interlace_t         interlace,
                                int                    n_parent_lists,
                                const cs_lnum_t        parent_num_shift[],
                                cs_datatype_t          datatype,
                                int                    time_step,
                                double                 time_value,
                                const void      *const field_values[])
{
  int  mesh_id, field_id;

  fvm_to_medcoupling_t *this_writer = (fvm_to_medcoupling_t *)this_writer_p;

  TypeOfTimeDiscretization  td = (time_step < 0) ? NO_TIME : ONE_TIME;

  /* Initialization */
  /*----------------*/

  mesh_id = _get_medcoupling_mesh_id(this_writer, mesh->name);

  if (mesh_id < 0) {
    mesh_id = _add_medcoupling_mesh(this_writer, mesh);
    fvm_to_medcoupling_export_nodal(this_writer, mesh);
  }

  /* Get field id */

  field_id = _get_medcoupling_field_id(this_writer,
                                       name,
                                       mesh_id,
                                       dimension,
                                       location,
                                       td);

  if (field_id < 0)
    field_id = _add_medcoupling_field(this_writer,
                                      name,
                                      mesh_id,
                                      dimension,
                                      location,
                                      td);

  MEDCouplingFieldDouble  *f = this_writer->fields[field_id]->f;

  /* Per node variable */
  /*-------------------*/

  if (location == FVM_WRITER_PER_NODE) {

    _export_field_values_n(mesh,
                           dimension,
                           interlace,
                           n_parent_lists,
                           parent_num_shift,
                           datatype,
                           field_values,
                           f);
  }

  /* Per element variable */
  /*----------------------*/

  else if (location == FVM_WRITER_PER_ELEMENT) {

    _export_field_values_e(mesh,
                           dimension,
                           interlace,
                           n_parent_lists,
                           parent_num_shift,
                           datatype,
                           field_values,
                           f);

  } /* End for per element variable */

  /* Update field status */
  /*---------------------*/

  if (td != NO_TIME)
    f->setTime(time_value, time_step, -1);
  f->getArray()->declareAsNew();

}

/*----------------------------------------------------------------------------
 * Flush files associated with a given writer.
 *
 * parameters:
 *   this_writer_p    <-- pointer to associated writer
 *----------------------------------------------------------------------------*/

void
fvm_to_medcoupling_flush(void  *this_writer_p)
{
  fvm_to_medcoupling_t *w = (fvm_to_medcoupling_t *)this_writer_p;

#if 0 /* example loop on fields */
  for (int i = 0; i < w->n_fields; i++) {
    w->fields[i]->f->...;
  }
#endif
}

/*----------------------------------------------------------------------------*/

} /* extern C */

#endif /* HAVE_MEDCOUPLING */
