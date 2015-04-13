/*============================================================================
 * Write a nodal representation associated with a mesh and associated
 * variables to CGNS files
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2013 EDF S.A.

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

#if defined(HAVE_CGNS)

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * CGNS library headers
 *----------------------------------------------------------------------------*/

#include <cgnslib.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_error.h>
#include <bft_mem.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_defs.h"
#include "fvm_convert_array.h"
#include "fvm_gather.h"
#include "fvm_io_num.h"
#include "fvm_nodal.h"
#include "fvm_nodal_priv.h"
#include "fvm_writer_helper.h"
#include "fvm_writer_priv.h"

#include "cs_parall.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvm_to_cgns.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

#define FVM_CGNS_NAME_SIZE      32      /* Maximum CGNS name length (the CGNS
                                           documentation does not specify this,
                                           but CGNS examples and ADF.h seem
                                           to use 32 character strings) */

/* Compatibility with different CGNS library versions */

#if !defined(CGNS_ENUMV)
#define CGNS_ENUMV(e) e
#endif

#if !defined(CGNS_ENUMT)
#define CGNS_ENUMT(e) e
#endif

#if !defined (CGNS_ENUMD)
#define CGNS_ENUMD(e) e
#endif

#if !defined (CGNS_ENUMF)
#define CGNS_ENUMF(e) e
#endif

#if CGNS_VERSION < 3100
#define cgsize_t int
#endif

/*============================================================================
 * Local Type Definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * CGNS solution structure
 *----------------------------------------------------------------------------*/

typedef struct {

  char                        *name;       /* Solution name */
  int                          index;      /* CGNS base index */
  CGNS_ENUMT(GridLocation_t)   location;   /* CGNS grid location of values */
  double                       time_value; /* Time step value */
  int                          time_step;  /* No. of iteration associated
                                              with time value */

} fvm_to_cgns_solution_t;

/*----------------------------------------------------------------------------
 * CGNS base structure
 *----------------------------------------------------------------------------*/

typedef struct {

  char     *name;              /* Mesh name */
  int       index;             /* CGNS base index */
  int       celldim;           /* Cell dimension:
                                  3 for a volume mesh, 2 for a surface mesh */
  int       physdim;           /* Physical dimension: number of coordinates
                                  defining a vertex */

  int       n_sols;                    /* Number of solutions */
  fvm_to_cgns_solution_t  **solutions; /* Array of pointers to CGNS solution
                                          structures in FVM */

} fvm_to_cgns_base_t;

/*----------------------------------------------------------------------------
 * CGNS writer structure
 *----------------------------------------------------------------------------*/

typedef struct {

  char                   *name;            /* Writer name */
  char                   *filename;        /* associated CGNS file name */
  int                     index;           /* index in associated CGNS file */

  int                     n_bases;         /* Number of CGNS bases */
  fvm_to_cgns_base_t    **bases;           /* Array of pointers to CGNS base
                                              structures in FVM */

  fvm_writer_time_dep_t   time_dependency; /* Mesh time dependency */
  int                     n_time_steps;    /* Number of mesh time steps */
  int                    *time_steps;      /* Array of mesh time steps */
  double                 *time_values;     /* Array of mesh time values */

  _Bool        is_open;            /* True if CGNS file is open */

  _Bool        discard_polygons;   /* Option to discard polygonal elements */
  _Bool        discard_polyhedra;  /* Option to discard polyhedral elements */

  _Bool        divide_polygons;    /* Option to tesselate polygonal elements */

  int          rank;               /* Rank of current process in communicator */
  int          n_ranks;            /* Number of processes in communicator */

#if defined(HAVE_MPI)
  MPI_Comm     comm;               /* Associated MPI communicator */
#endif

} fvm_to_cgns_writer_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

static char _cgns_version_string[32] = "";

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Delete CGNS base structure included in CGNS writer structure.
 *
 * parameters:
 *   base    <-- CGNS base structure in FVM.
 *----------------------------------------------------------------------------*/

static fvm_to_cgns_base_t *
_del_base(fvm_to_cgns_base_t  *base)
{
  int i;

  BFT_FREE(base->name);

  for (i = 0; i < base->n_sols; i++) {
    BFT_FREE(base->solutions[i]->name);
    BFT_FREE(base->solutions[i]);
  }

  BFT_FREE(base->solutions);
  BFT_FREE(base);

  return NULL;
}

/*----------------------------------------------------------------------------
 * Return the CGNS base index associated with a given CGNS base name, or 0.
 *
 * parameters:
 *   writer     <-- CGNS writer structure
 *   base_name  <-- base name
 *
 * returns:
 *    CGNS base index, or 0 if base name is not associated with this
 *    CGNS writer structure in FVM
 *----------------------------------------------------------------------------*/

static int
_get_base_index(fvm_to_cgns_writer_t  *writer,
                const char            *base_name)
{
  int i;

  fvm_to_cgns_base_t  **base_array = writer->bases;
  int retval = CG_OK;

  assert(writer != NULL);

  for (i = 0; i < writer->n_bases; i++) {
    if (strcmp(base_name, base_array[i]->name) == 0)
      break;
  }

  if (i == writer->n_bases)
    retval = 0;
  else
    retval = base_array[i]->index;

  return retval;
}

/*----------------------------------------------------------------------------
 * Associate a CGNS base name with a CGNS writer and return its index.
 * If the CGNS base was already associated, zero is returned.
 *
 * parameters:
 *   writer      <-- CGNS writer structure.
 *   base_name   <-- CGNS base name.
 *   mesh        <-- FVM mesh structure.
 *
 * returns:
 *   CGNS base index, or 0 if CGNS base already associated
 *----------------------------------------------------------------------------*/

static int
_add_base(fvm_to_cgns_writer_t  *writer,
          const char            *base_name,
          const fvm_nodal_t     *mesh)
{
  int  i;

  int  base_index = 0;
  int  rank = writer->rank;

  int entity_dim = fvm_nodal_get_max_entity_dim(mesh);

  int  retval = CG_OK;

  assert(writer != NULL);

  if (entity_dim == 0) /* CGNS requires entity_dim > 0 */
    entity_dim = mesh->dim;

  /* Add a new CGNS base structure */

  writer->n_bases += 1;
  i = writer->n_bases - 1;

  BFT_REALLOC(writer->bases, writer->n_bases, fvm_to_cgns_base_t *);
  BFT_MALLOC(writer->bases[i], 1, fvm_to_cgns_base_t);
  BFT_MALLOC(writer->bases[i]->name, strlen(base_name) + 1, char);

  strcpy(writer->bases[i]->name, base_name);
  writer->bases[i]->celldim  = entity_dim;
  writer->bases[i]->physdim  = mesh->dim;

  writer->bases[i]->n_sols = 0;
  writer->bases[i]->solutions = NULL;

  if (rank == 0) {
    retval = cg_base_write(writer->index,
                           base_name,
                           entity_dim,
                           mesh->dim,
                           &base_index);
    if (retval != CG_OK )
      bft_error(__FILE__, __LINE__, 0,
                _("cg_base_write() failed to create a new base:\n"
                  "Associated writer: \"%s\"\n"
                  "Associated mesh: \"%s\"\n%s"),
                writer->name, base_name, cg_get_error());
  }

#if defined(HAVE_MPI)
  if (writer->n_ranks > 1)
    MPI_Bcast(&base_index, 1, MPI_INT, 0, writer->comm);
#endif

  writer->bases[i]->index = base_index;

  return base_index;
}

/*----------------------------------------------------------------------------
 * Get base id associated with given index.
 *
 * parameters:
 *   writer      <-- pointer associated with writer
 *   base_index  <-- index of associated base
 *
 * returns:
 *   local base id for the given CGNS index
 *----------------------------------------------------------------------------*/

inline static int
_base_id(const fvm_to_cgns_writer_t   *writer,
         int                           base_index)
{
  int  base_id;

  assert(writer != NULL);
  assert(writer->bases != NULL);

  for (base_id = 0; base_id < writer->n_bases; base_id++) {
    fvm_to_cgns_base_t *base = writer->bases[base_id];
    if (base->index == base_index)
      return base_id;
  }

  bft_error(__FILE__, __LINE__, 0,
            _("No CGNS base with index %d defined:\n"
              "Associated writer: \"%s\"\n"),
            base_index, writer->name);
  return -1;
}

/*----------------------------------------------------------------------------
 * Associate new time step with a CGNS writer.
 *
 * parameters:
 *   writer      <-- pointer associated with writer
 *   base_index  <-- index of associated base
 *   time_step   <-- time step number
 *   time_value  <-- time_value number
 *   location    <-- CGNS grid location
 *
 * returns:
 *   solution index of the new time step
 *----------------------------------------------------------------------------*/

static int
_add_solution(fvm_to_cgns_writer_t        *writer,
              int                          base_index,
              int                          time_step,
              double                       time_value,
              CGNS_ENUMT(GridLocation_t)   location)
{
  int  sol_id, sol_length;
  char sol_name[FVM_CGNS_NAME_SIZE + 1];

  int  zone_index = 1; /* We always write to the first zone */
  int  sol_index = 0;
  int  rank = writer->rank;

  int  retval = CG_OK;

  /* Find matching base */

  int base_id = _base_id(writer, base_index);

  /* Create a new pointer to fvm_to_cgns_solution_t */

  fvm_to_cgns_base_t *base = writer->bases[base_id];

  base->n_sols += 1;
  sol_id = base->n_sols - 1;

  BFT_REALLOC(base->solutions, sol_id + 1, fvm_to_cgns_solution_t *);
  BFT_MALLOC(base->solutions[sol_id], 1, fvm_to_cgns_solution_t);

  /* Initialization of the new solution structure */

  base->solutions[sol_id]->index = -1;
  base->solutions[sol_id]->time_step = time_step;
  base->solutions[sol_id]->time_value = time_value;
  base->solutions[sol_id]->location = location;
  base->solutions[sol_id]->name = NULL;

  if (time_step < 0)
    sprintf(sol_name, "Stationary (%s)",
            GridLocationName[location]);
  else
    sprintf(sol_name, "Solution %3d (%s)",
            (int)time_step, GridLocationName[location]);

  sol_length = strlen(sol_name);
  BFT_MALLOC(base->solutions[sol_id]->name, sol_length + 1, char);
  strncpy(base->solutions[sol_id]->name, sol_name, sol_length);
  base->solutions[sol_id]->name[sol_length] = '\0';

  if (rank == 0) {
    retval = cg_sol_write(writer->index,
                          base->index,
                          zone_index,
                          base->solutions[sol_id]->name,
                          location,
                          &sol_index);

    if (retval != CG_OK)
      bft_error(__FILE__, __LINE__, 0,
                _("cg_sol_write() failed to create a "
                  "new solution node:\n"
                  "Associated writer: \"%s\"\n"
                  "Associated mesh: \"%s\"\n"
                  "Solution name: \"%s\"\n"
                  "Associated time value: %f \n%s"),
                writer->name, base->name,
                base->solutions[sol_id]->name, time_value, cg_get_error());
  }

#if defined(HAVE_MPI)
  if (writer->n_ranks > 1)
    MPI_Bcast(&sol_index, 1, MPI_INT, 0, writer->comm);
#endif

  base->solutions[sol_id]->index = sol_index;

  return (sol_index);
}

/*----------------------------------------------------------------------------
 * Find the solution index associated with a given time step.
 *
 * parameters:
 *   writer      <-- pointer to associated writer
 *   base_index  <-- index of associated base
 *   time_step   <-- time step number
 *   time_value  <-- time_value number
 *   location    <-- location of results (CellCenter, Vertex, ...)
 *
 * returns:
 *   solution index associated with given time step, or 0 if none found.
 *----------------------------------------------------------------------------*/

static int
_get_solution_index(const fvm_to_cgns_writer_t  *writer,
                    int                          base_index,
                    int                          time_step,
                    double                       time_value,
                    CGNS_ENUMT(GridLocation_t)   location)
{
  int sol_id, n_sols;

  int sol_index = 0;
  fvm_to_cgns_base_t     **bases = writer->bases;
  fvm_to_cgns_solution_t  *sol_ref = NULL;

  const char time_value_err_string[] =
    N_("The time value associated with time step <%d> equals <%g>,\n"
       "but time value <%g> has already been associated with this time step.\n");

  /* Find matching base */

  int base_id = _base_id(writer, base_index);

  /* Any negative time step value indicates time independant values */

  if (time_step < 0) {
    time_step = -1;
    time_value = 0.0;
  }

  if (bases[base_id]->solutions != NULL) {
    n_sols = bases[base_id]->n_sols;

    /* Search for index associated with time step */

    for (sol_id = 0; sol_id < n_sols; sol_id++) {

      sol_ref = bases[base_id]->solutions[sol_id];

      /* Check on grid location */
      if (location == sol_ref->location) {

        /* Check on time step */
        if (time_step == sol_ref->time_step) {

        /* Check on time value */
          if (time_value < sol_ref->time_value - 1.e-16 ||
              time_value > sol_ref->time_value + 1.e-16)
            bft_error(__FILE__, __LINE__, 0,
                _(time_value_err_string), time_step,
                time_value, sol_ref->time_value);

          else {
            sol_index = sol_ref->index;
            break;
          }
        }
        else
          sol_index = 0;
      }
      else
        sol_index = 0;

    } /* End of loop on existing solutions */

  }
  else  /* Set the first solution */
    sol_index = 0;

  return sol_index;
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
_count_extra_vertices(const fvm_to_cgns_writer_t *this_writer,
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

    const fvm_nodal_section_t  *const  section = mesh->sections[i];

    /* Output if entity dimension equal to highest in mesh
       (i.e. no output of faces if cells present, or edges
       if cells or faces) */

    if (   section->entity_dim == export_dim
        && section->type == FVM_CELL_POLY
        && section->tesselation != NULL
        && this_writer->discard_polyhedra == false) {

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
 * Create a CGNS zone associated with the writer.
 * Only one zone is created per base. Its index must be 1.
 *
 * parameters:
 *   mesh                 <-- FVM mesh structure.
 *   writer               <-- CGNS writer structure.
 *   base                 <-- CGNS base structure.
 *   export_sections      <-> pointer to a list of fvm_writer_section structures
 *   n_g_vertices         <-- global number of vertices in the mesh
 *
 *----------------------------------------------------------------------------*/

static void
_add_zone(const fvm_nodal_t           *mesh,
          const fvm_to_cgns_writer_t  *writer,
          const fvm_to_cgns_base_t    *base,
          const fvm_writer_section_t  *export_sections,
          cs_gnum_t                    n_g_vertices)
{
  int zone_index;
  cgsize_t    zone_sizes[3];

  cs_gnum_t   n_g_entities = 0;
  cs_gnum_t   n_g_tesselated_elements = 0;
  cs_gnum_t   n_g_extra_vertices = 0;

  const fvm_writer_section_t *current_section = NULL;

  int retval = CG_OK;

  assert(writer != NULL);

  if (writer->rank != 0)
    return;

  /* Compute global number of vertices in this zone */

  /* Polyhedra are not defined in CGNS. If they are not discarded,
     they have to be tesselated. */

  _count_extra_vertices(writer,
                        mesh,
                        &n_g_extra_vertices,
                        NULL);


  zone_sizes[0] = n_g_vertices + n_g_extra_vertices;

  /* Compute global number of entities in this zone */

  current_section = export_sections;

  while (current_section != NULL) {

    const fvm_nodal_section_t *const section = current_section->section;

    if (current_section->type == section->type)

      /* Regular section */
      n_g_entities += fvm_nodal_section_n_g_elements(section);

    else {

      /* Tesselated section */
      fvm_tesselation_get_global_size(section->tesselation,
                                      current_section->type,
                                      &n_g_tesselated_elements,
                                      NULL);

      n_g_entities += n_g_tesselated_elements;
    }

    current_section = current_section->next;

  }  /* End of loop on sections */

  zone_sizes[1] = n_g_entities;

  /* Set boundary vertex size (zero if element not sorted) */

  zone_sizes[2] = 0;

  /* Create CGNS zone */

  assert(writer->is_open == true);

  retval = cg_zone_write(writer->index,
                         base->index,
                         "Zone 1",
                         zone_sizes,
                         CGNS_ENUMV(Unstructured),
                         &zone_index);

  if (retval != CG_OK)
    bft_error(__FILE__, __LINE__, 0,
              _("cg_zone_write() failed to create a new zone:\n"
                "Associated writer: \"%s\"\n"
                "Associated base: \"%s\"\n%s"),
              writer->name, base->name, cg_get_error());

  assert(zone_index == 1);

  return;
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
_extra_vertex_coords(const fvm_to_cgns_writer_t  *this_writer,
                     const fvm_nodal_t           *mesh)
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

      const fvm_nodal_section_t  *const  section = mesh->sections[i];

      /* Output if entity dimension equal to highest in mesh
         (i.e. no output of faces if cells present, or edges
         if cells or faces) */

      if (   section->entity_dim == export_dim
          && section->type == FVM_CELL_POLY
          && section->tesselation != NULL
          && this_writer->discard_polyhedra == false) {

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
 * Define section name and CGNS element type.
 *
 * parameters:
 *   export_section  <-- pointer to section list structure.
 *   section_id      <-- identification number of the section.
 *   section_name    --> name of section for CGNS.
 *   cgns_elem_type  --> CGNS element type for this section.
 *----------------------------------------------------------------------------*/

static void
_define_section(const fvm_writer_section_t  *section,
                int section_id,
                char section_name[FVM_CGNS_NAME_SIZE + 1],
                CGNS_ENUMT(ElementType_t) *cgns_elt_type)
{
  assert(section != NULL);

  switch(section->type) {

  case FVM_EDGE:
    sprintf(section_name, "Edges_%d", section_id);
    *cgns_elt_type = CGNS_ENUMV(BAR_2);
    break;

  case FVM_FACE_TRIA:
    sprintf(section_name, "Triangles_%d", section_id);
    *cgns_elt_type = CGNS_ENUMV(TRI_3);
    break;

  case FVM_FACE_QUAD:
    sprintf(section_name, "Quadrangles_%d", section_id);
    *cgns_elt_type = CGNS_ENUMV(QUAD_4);
    break;

  case FVM_FACE_POLY:
    sprintf(section_name, "Polygons_%d", section_id);
#if CGNS_VERSION < 3200
    *cgns_elt_type = CGNS_ENUMV(MIXED);
#else
    *cgns_elt_type = CGNS_ENUMV(NGON_n);
#endif
    break;

  case FVM_CELL_TETRA:
    sprintf(section_name, "Tetrahedra_%d", section_id);
    *cgns_elt_type = CGNS_ENUMV(TETRA_4);
    break;

  case FVM_CELL_PYRAM:
    sprintf(section_name, "Pyramids_%d", section_id);
    *cgns_elt_type = CGNS_ENUMV(PYRA_5);
    break;

  case FVM_CELL_PRISM:
    sprintf(section_name, "Prisms_%d", section_id);
    *cgns_elt_type = CGNS_ENUMV(PENTA_6);
    break;

  case FVM_CELL_HEXA:
    sprintf(section_name, "Hexahedra_%d", section_id);
    *cgns_elt_type = CGNS_ENUMV(HEXA_8);
    break;

  default:
    sprintf(section_name, "Null_section_%d", section_id);
    *cgns_elt_type = CGNS_ENUMV(ElementTypeNull);
  }

  return;
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Write vertex coordinates to a CGNS file in parallel mode
 *
 * parameters:
 *   writer        <-- pointer to associated writer.
 *   mesh          <-- pointer to nodal mesh structure that should be written.
 *   base          <-- pointer to CGNS base structure.
 *   global_s_size <-- maximum number of entities defined per slice.
 *----------------------------------------------------------------------------*/

static void
_export_vertex_coords_g(fvm_to_cgns_writer_t  *writer,
                        const fvm_nodal_t     *mesh,
                        fvm_to_cgns_base_t    *base,
                        cs_gnum_t              global_s_size)
{
  int coord_index, section_id;
  cs_lnum_t   i, j;
  MPI_Datatype  mpi_datatype;
  CGNS_ENUMT(DataType_t)  cgns_datatype;

  cgsize_t  idx_start = 0, idx_end = 0;
  int zone_index = 1;
  cs_lnum_t   n_extra_vertices = 0;
  cs_gnum_t   n_g_extra_vertices = 0;
  cs_gnum_t   global_num_start = 0, global_num_end = 0;
  cs_coord_t *extra_vertex_coords = NULL;
  cs_coord_t *coords_tmp = NULL;
  cs_coord_t *global_coords_s = NULL;
  fvm_gather_slice_t  *vertices_slice = NULL;

  const cs_coord_t *vertex_coords = mesh->vertex_coords;
  const cs_lnum_t *parent_vertex_num = mesh->parent_vertex_num;
  const cs_lnum_t n_vertices = mesh->n_vertices;
  const cs_gnum_t n_g_vertices
    = fvm_io_num_get_global_count(mesh->global_vertex_num);

  const char *const coord_name[3] = {"CoordinateX",
                                     "CoordinateY",
                                     "CoordinateZ"};

  size_t  stride = (size_t)mesh->dim;
  int  retval = CG_OK;

  assert(base != NULL);

  if (sizeof(cs_coord_t) == sizeof(double)) {
    cgns_datatype = CGNS_ENUMV(RealDouble);
    mpi_datatype = MPI_DOUBLE;
  }
  else {
    cgns_datatype = CGNS_ENUMV(RealSingle);
    mpi_datatype = MPI_FLOAT;
  }

  _count_extra_vertices(writer,
                        mesh,
                        &n_g_extra_vertices,
                        &n_extra_vertices);

  extra_vertex_coords = _extra_vertex_coords(writer,
                                             mesh);

  BFT_MALLOC(coords_tmp, CS_MAX(n_vertices, n_extra_vertices), cs_coord_t);
  BFT_MALLOC(global_coords_s, global_s_size, cs_coord_t);

  vertices_slice = fvm_gather_slice_create(mesh->global_vertex_num,
                                           global_s_size,
                                           writer->comm);

  /* Loop on spatial dimension */

  for (j = 0; j < base->physdim; j ++) {

    /* Vertex coordinates */

    if (parent_vertex_num != NULL) {
      for (i = 0; i < n_vertices; i++)
        coords_tmp[i] = vertex_coords[(parent_vertex_num[i]-1)*stride + j];
    }
    else {
      for (i = 0; i < n_vertices; i++)
        coords_tmp[i] = vertex_coords[i*stride + j];
    }

    if (j > 0)
      fvm_gather_slice_reinitialize(vertices_slice);

    while (fvm_gather_slice_advance(vertices_slice,
                                    &global_num_start,
                                    &global_num_end) == 0) {

      fvm_gather_array(coords_tmp,
                       global_coords_s,
                       mpi_datatype,
                       1,
                       mesh->global_vertex_num,
                       writer->comm,
                       vertices_slice);

      if (writer->rank == 0) {        /* Write grid coordinates */

        if (global_coords_s != NULL) {

          idx_start = global_num_start;
          idx_end = global_num_end - 1;

          retval = cg_coord_partial_write(writer->index,
                                          base->index,
                                          zone_index,
                                          cgns_datatype,
                                          coord_name[j],
                                          &idx_start,
                                          &idx_end,
                                          global_coords_s,
                                          &coord_index);

        }

        if (retval != CG_OK)
          bft_error(__FILE__, __LINE__, 0,
                    _("cg_coord_partial_write() failed to write coords:\n"
                      "Associated writer: \"%s\"\n"
                      "Associated base: \"%s\"\n"
                      "CGNS error:%s"),
                    writer->name, base->name, cg_get_error());

      } /* End if rank == 0 */

    } /* End of slice advance */

    /* Extra vertex coordinates */

    if (n_g_extra_vertices > 0) {

      cs_lnum_t extra_vertices_count = 0;

      for (section_id = 0 ; section_id < mesh->n_sections ; section_id++) {

        fvm_gather_slice_t  *extra_vertices_slice = NULL;
        const fvm_nodal_section_t  *const  section = mesh->sections[section_id];

        /* Output if entity dimension equal to highest in mesh
           (i.e. no output of faces if cells present, or edges
           if cells or faces) */

        if (   section->entity_dim == mesh->dim
            && section->type == FVM_CELL_POLY
            && section->tesselation != NULL
            && writer->discard_polyhedra == false) {

          const fvm_io_num_t *extra_vertex_num
            = fvm_tesselation_global_vertex_num(section->tesselation);
          const cs_lnum_t n_extra_vertices_section
            = fvm_tesselation_n_vertices_add(section->tesselation);

          for (i = 0 ; i < n_extra_vertices ; i++)
            coords_tmp[i]
              = (extra_vertex_coords[(i+extra_vertices_count)*stride + j]);

          extra_vertices_slice = fvm_gather_slice_create(extra_vertex_num,
                                                         global_s_size,
                                                         writer->comm);

          /* loop on slices in parallel mode */

          while (fvm_gather_slice_advance(extra_vertices_slice,
                                          &global_num_start,
                                          &global_num_end) == 0) {

            fvm_gather_array(coords_tmp,
                             global_coords_s,
                             mpi_datatype,
                             1,
                             extra_vertex_num,
                             writer->comm,
                             extra_vertices_slice);

            if (writer->rank == 0) {

              if (global_coords_s != NULL) {

                idx_start = global_num_start + n_g_vertices;
                idx_end = global_num_end - 1 + n_g_vertices;

                retval = cg_coord_partial_write(writer->index,
                                                base->index,
                                                zone_index,
                                                cgns_datatype,
                                                coord_name[j],
                                                &idx_start,
                                                &idx_end,
                                                global_coords_s,
                                                &coord_index);

              }

              if (retval != CG_OK)
                bft_error(__FILE__, __LINE__, 0,
                          _("cg_coord_partial_write() failed to write "
                            "extra vertices coords:\n"
                            "Associated writer: \"%s\"\n"
                            "Associated base: \"%s\"\n"
                            "CGNS error:%s"),
                          writer->name, base->name, cg_get_error());

            } /* End if rank == 0 */

          } /* End of loop on slices */

          fvm_gather_slice_destroy(extra_vertices_slice);

          extra_vertices_count += n_extra_vertices_section;

        }

      } /* end of loop on sections for extra vertices */

    } /* End of treatment of extra vertices */

  } /* End of loop on spatial dimension */

  fvm_gather_slice_destroy(vertices_slice);

  BFT_FREE(coords_tmp);
  BFT_FREE(global_coords_s);

  if (extra_vertex_coords != NULL)
    BFT_FREE(extra_vertex_coords);

  return;
}

#endif

/*----------------------------------------------------------------------------
 * Write vertex coordinates to a CGNS file in serial mode
 *
 * parameters:
 *   writer        <-- pointer to associated writer.
 *   mesh          <-- pointer to nodal mesh structure.
 *   base          <-- pointer to CGNS base structure.
 *----------------------------------------------------------------------------*/

static void
_export_vertex_coords_l(const fvm_to_cgns_writer_t  *writer,
                        const fvm_nodal_t     *mesh,
                        fvm_to_cgns_base_t    *base)
{
  int  coord_index;
  cs_lnum_t   i, j;
  cgsize_t  idx_start, idx_end;
  CGNS_ENUMT(DataType_t)  cgns_datatype;

  int  zone_index = 1;
  size_t  stride = (size_t)mesh->dim;
  cs_lnum_t   n_extra_vertices = 0;
  cs_coord_t  *extra_vertex_coords = NULL;
  cs_coord_t  *coords_tmp = NULL;

  const cs_lnum_t   n_vertices = mesh->n_vertices;
  const cs_lnum_t *parent_vertex_num = mesh->parent_vertex_num;
  const cs_coord_t *vertex_coords = mesh->vertex_coords;

  const char *const coord_name[3] = {"CoordinateX",
                                     "CoordinateY",
                                     "CoordinateZ"};

  int  retval = CG_OK;

  assert(writer->is_open == true);
  assert(base != NULL);

  if (sizeof(cs_coord_t) == sizeof(double))
    cgns_datatype = CGNS_ENUMV(RealDouble);
  else
    cgns_datatype = CGNS_ENUMV(RealSingle);

  /* Compute extra vertex coordinates if present */

  _count_extra_vertices(writer,
                        mesh,
                        NULL,
                        &n_extra_vertices);

  extra_vertex_coords = _extra_vertex_coords(writer,
                                             mesh);


  BFT_MALLOC(coords_tmp, CS_MAX(n_vertices, n_extra_vertices), cs_coord_t);

  /* Loop on dimension */

  for (j = 0; j < base->physdim; j ++) {

    /* Vertex coordinates */

    if (parent_vertex_num != NULL) {
      for (i = 0; i < n_vertices; i++)
        coords_tmp[i] = vertex_coords[(parent_vertex_num[i]-1)*stride + j];
    }
    else {
      for (i = 0; i < n_vertices; i++)
        coords_tmp[i] = vertex_coords[i*stride + j];
    }

    /* Write grid coordinates */

    if (coords_tmp != NULL) {

      idx_start = 1;
      idx_end = mesh->n_vertices;

      retval = cg_coord_partial_write(writer->index,
                                      base->index,
                                      zone_index,
                                      cgns_datatype,
                                      coord_name[j],
                                      &idx_start,
                                      &idx_end,
                                      coords_tmp,
                                      &coord_index);

    }

    if (retval != CG_OK)
      bft_error(__FILE__, __LINE__, 0,
                _("cg_coord_partial write() failed to write coords:\n"
                  "Associated writer: \"%s\"\n"
                  "Associated base: \"%s\"\n"
                  "Associated zone index: \"%i\"\n%s"),
                writer->name, base->name, zone_index, cg_get_error());

    /* Extra vertex coordinates */

    for (i = 0 ; i < n_extra_vertices ; i++)
      coords_tmp[i] = extra_vertex_coords[i*stride + j];

    if (n_extra_vertices > 0) {

      idx_start = mesh->n_vertices + 1;
      idx_end = mesh->n_vertices + n_extra_vertices;

      retval = cg_coord_partial_write(writer->index,
                                      base->index,
                                      zone_index,
                                      cgns_datatype,
                                      coord_name[j],
                                      &idx_start,
                                      &idx_end,
                                      coords_tmp,
                                      &coord_index);

      if (retval != CG_OK)
        bft_error(__FILE__, __LINE__, 0,
                  _("cg_coord_partial_write() failed to write"
                    " extra vertices coords:\n"
                    "Associated writer: \"%s\"\n"
                    "Associated base: \"%s\"\n"
                    "Associated zone index: \"%i\"\n%s"),
                  writer->name, base->name, zone_index, cg_get_error());

    }

  } /* End of loop on coordinates */

  BFT_FREE(coords_tmp);

  if (extra_vertex_coords != NULL)
    BFT_FREE(extra_vertex_coords);

}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Write strided connectivity to a CGNS file in parallel mode
 *
 * parameters:
 *   export_section   <-- pointer to sections list to export.
 *   writer           <-- pointer to associated writer.
 *   mesh             <-- pointer to nodal mesh structure.
 *   base             <-- pointer to CGNS base structure.
 *   section_id       <-- section identificator number.
 *   global_counter   <-- counter to update element shift after each section.
 *   global_s_size    <-- maximum number of entities defined per slice
 *   global_connect_s  --> global connectivity array section for elements
 *                         slice global_num_start to global_num_end
 *                         (output for rank 0, working array only for others)
 *----------------------------------------------------------------------------*/

static const fvm_writer_section_t *
_export_connect_g(const fvm_writer_section_t  *export_section,
                  fvm_to_cgns_writer_t        *writer,
                  const fvm_nodal_t           *mesh,
                  fvm_to_cgns_base_t          *base,
                  int                          section_id,
                  cs_gnum_t                   *global_counter,
                  cs_gnum_t                    global_s_size,
                  cs_gnum_t                    global_connect_s_call[])
{
  char  section_name[FVM_CGNS_NAME_SIZE + 1];
  CGNS_ENUMT(ElementType_t) cgns_elt_type; /* Definition in cgnslib.h */

  int  section_index = -1;
  cs_gnum_t   global_num_start = 0, global_num_end = 0;
  cs_gnum_t   elt_start = 0, elt_end = 0;
  cs_gnum_t   *global_connect_s = global_connect_s_call;
  fvm_gather_slice_t  *elements_slice = NULL;

  const int  zone_index = 1; /* We always use zone index = 1 */
  const fvm_writer_section_t *current_section = export_section;
  const fvm_nodal_section_t *section = current_section->section;
  const fvm_io_num_t  *global_element_num = section->global_element_num;

  int  retval = CG_OK;

  assert(current_section != NULL);

  _define_section(current_section, section_id, section_name, &cgns_elt_type);

  elements_slice
    = fvm_gather_slice_create(global_element_num,
                              global_s_size,
                              writer->comm);

  while (fvm_gather_slice_advance(elements_slice,
                                  &global_num_start,
                                  &global_num_end) == 0) {

    fvm_gather_strided_connect(section->vertex_num,
                               global_connect_s,
                               section->stride,
                               mesh->global_vertex_num,
                               global_element_num,
                               writer->comm,
                               elements_slice);

    if (writer->rank == 0) {

      elt_start = *global_counter + global_num_start;
      elt_end   = *global_counter + global_num_end - 1;

      if (global_connect_s != NULL) {
        cgsize_t *_global_connect_s = (cgsize_t *)global_connect_s;
        if (sizeof(cgsize_t) != sizeof(cs_gnum_t)) {
          cgsize_t i = 0, n = (elt_end + 1 - elt_start)*section->stride;
          if (sizeof(cgsize_t) > sizeof(cs_gnum_t))
            BFT_MALLOC(_global_connect_s, n, cgsize_t);
          for (i = 0; i < n; i++)
            _global_connect_s[i] = global_connect_s[i];
        }

        if (global_num_start == 1) { /* First pass */
          retval = cg_section_partial_write(writer->index,
                                            base->index,
                                            zone_index,
                                            section_name,
                                            cgns_elt_type,
                                            elt_start,
                                            elt_end,
                                            0, /* unsorted boundary elements */
                                            &section_index);
          if (retval != CG_OK)
            bft_error(__FILE__, __LINE__, 0,
                      _("cg_section_partial_write() failed to write elements:\n"
                        "Associated writer: \"%s\"\n"
                        "Associated base: \"%s\"\n"
                        "Associated section name: \"%s\"\n%s"),
                      writer->name, base->name, section_name, cg_get_error());
        }

        if (retval == CG_OK)
          retval = cg_elements_partial_write(writer->index,
                                             base->index,
                                             zone_index,
                                             section_index,
                                             elt_start,
                                             elt_end,
                                             _global_connect_s);
        if (retval != CG_OK)
          bft_error(__FILE__, __LINE__, 0,
                    _("cg_elements_partial_write() failed to write elements:\n"
                      "Associated writer: \"%s\"\n"
                      "Associated base: \"%s\"\n"
                      "Associated section name: \"%s\"\n"
                      "Associated range: [%llu, %llu]\n%s\n"),
                    writer->name, base->name, section_name,
                    (unsigned long long) elt_start,
                    (unsigned long long) elt_end,
                    cg_get_error());

        if (sizeof(cgsize_t) > sizeof(cs_gnum_t))
          BFT_FREE(_global_connect_s);
      }


    } /* End if rank = 0 */

  } /* End of loop on slice advance */

  if (elt_end >= elt_start)
    *global_counter += global_num_end - 1;

  /* Free slice structure */

  fvm_gather_slice_destroy(elements_slice);

  current_section = current_section->next;

  return current_section;
}
#endif

/*----------------------------------------------------------------------------
 * Write strided connectivity to a CGNS file in serial mode
 *
 * parameters:
 *   export_section   <-- pointer to sections list to export.
 *   writer           <-- pointer to associated writer.
 *   base             <-- pointer to CGNS base structure.
 *   section_id       <-- section identificator number.
 *   global_counter   <-- counter to update the shift after each section export.
 *----------------------------------------------------------------------------*/

static const fvm_writer_section_t *
_export_connect_l(const fvm_writer_section_t  *export_section,
                  const fvm_to_cgns_writer_t  *writer,
                  const fvm_to_cgns_base_t    *base,
                  int                          section_id,
                  cs_gnum_t                   *global_counter)
{
  int  section_index;
  char section_name[FVM_CGNS_NAME_SIZE + 1];
  CGNS_ENUMT(ElementType_t) cgns_elt_type; /* Definition in cgnslib.h */

  cs_gnum_t elt_start = 0, elt_end = 0;

  const int  zone_index = 1; /* We always use zone index = 1 */
  const fvm_writer_section_t  *current_section = export_section;
  const fvm_nodal_section_t *const section = current_section->section;

  int  retval = CG_OK;

  assert(current_section != NULL);

  _define_section(current_section, section_id, section_name, &cgns_elt_type);

  elt_start = 1 + *global_counter;
  elt_end   = section->n_elements + *global_counter;

  if (section->vertex_num != NULL) {
    cgsize_t *_vertex_num = NULL;
    const cgsize_t *vertex_num = (const cgsize_t *)section->vertex_num;
    if (sizeof(cgsize_t) != sizeof(cs_lnum_t)) {
      cgsize_t i = 0, n = (elt_end + 1 - elt_start)*section->stride;
      BFT_MALLOC(_vertex_num, n, cgsize_t);
      vertex_num = _vertex_num;
      for (i = 0; i < n; i++)
        _vertex_num[i] = section->vertex_num[i];
    }
    retval = cg_section_write(writer->index,
                              base->index,
                              zone_index,
                              section_name,
                              cgns_elt_type,
                              elt_start,
                              elt_end,
                              0, /* unsorted boundary elements */
                              vertex_num,
                              &section_index);
    if (_vertex_num != NULL)
      BFT_FREE(_vertex_num);
  }

  if (retval != CG_OK)
    bft_error(__FILE__, __LINE__, 0,
              _("cg_section_write() failed to write elements:\n"
                "Associated writer: \"%s\"\n"
                "Associated base: \"%s\"\n"
                "Associated section name: \"%s\"\n%s"),
              writer->name, base->name, section_name, cg_get_error());

  *global_counter += section->n_elements;

  current_section = current_section->next;

  return current_section;
}

#if defined(HAVE_MPI)
/*----------------------------------------------------------------------------
 * Write strided connectivity from tesselated elements to a CGNS file
 * in parallel mode.
 *
 * parameters:
 *   export_section   <-- pointer to sections list to export.
 *   writer           <-- pointer to associated writer.
 *   base             <-- pointer to CGNS base structure.
 *   section_id       <-- section identificator number.
 *   global_counter   <-- counter to update the shift after each section export.
 *   global_s_size    <-- maximum number of entities defined per slice.
 *   global_connect_s_size_call <-- buffer size to export connectivity
 *   global_connect_s <-- global connectivity array section for elements
 *                        slice global_num_start to global_num_end.
 *                        (output for rank 0, working array only for others)
 *----------------------------------------------------------------------------*/

static const fvm_writer_section_t *
_export_nodal_tesselated_g(const fvm_writer_section_t  *export_section,
                           const fvm_to_cgns_writer_t  *writer,
                           const fvm_nodal_t  *mesh,
                           const fvm_to_cgns_base_t  *base,
                           int  section_id,
                           cs_gnum_t   *global_counter,
                           cs_gnum_t   global_s_size,
                           cs_gnum_t   global_connect_s_size_call,
                           cs_gnum_t   global_connect_s_call[])
{
  char section_name[FVM_CGNS_NAME_SIZE + 1];
  CGNS_ENUMT(ElementType_t) cgns_elt_type; /* Definition in cgnslib.h */

  int  section_index = -1;
  cs_lnum_t   n_sub_elements_max = 0;
  cs_lnum_t   start_id = 0, end_id = 0;
  cs_gnum_t   elt_start = 0, elt_end = 0;
  cs_gnum_t   global_num_start = 0, global_num_end = 0;
  cs_gnum_t   n_g_sub_elements = 0, local_connect_size = 0;
  cs_gnum_t   global_connect_s_size_prev = global_connect_s_size_call;
  cs_gnum_t   global_connect_s_size = global_connect_s_size_call;

  cs_lnum_t   *local_idx = NULL;
  cs_gnum_t   *global_idx_s = NULL;
  cs_gnum_t   *sub_elt_vertex_num = NULL;
  cs_gnum_t   *global_connect_s = global_connect_s_call;

  fvm_gather_slice_t  *elements_slice = NULL;

  const int  zone_index = 1; /* We always use zone index = 1 */
  const fvm_writer_section_t *current_section = export_section;
  const fvm_nodal_section_t *section = current_section->section;
  const fvm_tesselation_t *tesselation = section->tesselation;
  const cs_gnum_t extra_vertex_base = current_section->extra_vertex_base;
  const cs_lnum_t n_elements = fvm_tesselation_n_elements(tesselation);
  const int stride = fvm_nodal_n_vertices_element[current_section->type];

  int  retval = CG_OK;

  assert(current_section != NULL);

  _define_section(current_section, section_id, section_name, &cgns_elt_type);

  fvm_tesselation_get_global_size(tesselation,
                                  current_section->type,
                                  &n_g_sub_elements,
                                  &n_sub_elements_max);

  /* Allocate memory for additionnal indexes and decoded connectivity */

  BFT_MALLOC(local_idx, n_elements + 1, cs_lnum_t);
  BFT_MALLOC(global_idx_s, global_s_size + 1, cs_gnum_t);

  local_connect_size = CS_MAX(global_s_size,
                              (cs_gnum_t)n_sub_elements_max*10);
  BFT_MALLOC(sub_elt_vertex_num, local_connect_size * stride, cs_gnum_t);

  /* Loop on slices */
  /*----------------*/

  /* fvm_tesselation_dump(tesselation); */

  elements_slice = fvm_gather_slice_create(section->global_element_num,
                                           local_connect_size,
                                           writer->comm);

  while (fvm_gather_slice_advance(elements_slice,
                                  &global_num_start,
                                  &global_num_end) == 0) {

    /* Build element->vertices index */

    end_id
      = fvm_tesselation_range_index_g(tesselation,
                                      current_section->type,
                                      stride,
                                      start_id,
                                      local_connect_size,
                                      &global_num_end,
                                      local_idx,
                                      writer->comm);

    /* Check if the maximum id returned on some ranks leads to a
       lower global_num_end than initially required (due to the
       local buffer being too small) and adjust slice if necessary */

    fvm_gather_slice_limit(elements_slice, &global_num_end);

    /* Gather element->vertices index */

    fvm_gather_slice_index(local_idx,
                           global_idx_s,
                           section->global_element_num,
                           writer->comm,
                           elements_slice);

    /* Recompute maximum value of global_num_end for this slice */

    fvm_gather_resize_indexed_slice(10,
                                    &global_num_end,
                                    &global_connect_s_size,
                                    writer->comm,
                                    global_idx_s,
                                    elements_slice);

    /* If the buffer passed to this function is too small, allocate a
       larger one; in this case, we may as well keep it for all slices */

    if (global_connect_s_size_prev < global_connect_s_size) {
      if (global_connect_s == global_connect_s_call)
        global_connect_s = NULL;
      BFT_REALLOC(global_connect_s, global_connect_s_size, cs_gnum_t);
      global_connect_s_size_prev = global_connect_s_size;
    }

    /* Now decode tesselation */

    end_id = fvm_tesselation_decode_g(tesselation,
                                      current_section->type,
                                      start_id,
                                      local_connect_size,
                                      &global_num_end,
                                      mesh->global_vertex_num,
                                      extra_vertex_base,
                                      sub_elt_vertex_num,
                                      writer->comm);

    /* No need to check if the maximum id returned on some ranks
       leads to a lower global_num_end than initially required
       (due to local buffer being full), as this was already done
       above for the local index */

    /* Now gather decoded element->vertices connectivity */

    fvm_gather_indexed(sub_elt_vertex_num,
                       global_connect_s,
                       CS_MPI_GNUM,
                       local_idx,
                       section->global_element_num,
                       writer->comm,
                       global_idx_s,
                       elements_slice);

    /* Do all printing for cells on rank 0 */

    if (writer->rank == 0) {

      elt_start = *global_counter + 1 + global_idx_s[0] / stride;
      elt_end   = *global_counter +
        global_idx_s[(global_num_end - global_num_start)] / stride;

      if (global_connect_s != NULL) {
        cgsize_t *_global_connect_s = (cgsize_t *)global_connect_s;
        if (sizeof(cgsize_t) != sizeof(cs_gnum_t)) {
          int i = 0, n = (elt_end + 1 - elt_start)*stride;
          if (sizeof(cgsize_t) > sizeof(cs_gnum_t))
            BFT_MALLOC(_global_connect_s, n, cgsize_t);
          for (i = 0; i < n; i++)
            _global_connect_s[i] = global_connect_s[i];
        }

        if (start_id == 0) { /* First pass */
          retval = cg_section_partial_write(writer->index,
                                            base->index,
                                            zone_index,
                                            section_name,
                                            cgns_elt_type,
                                            elt_start,
                                            elt_end,
                                            0, /* unsorted boundary elements */
                                            &section_index);
          if (retval != CG_OK)
            bft_error(__FILE__, __LINE__, 0,
                      _("cg_section_partial_write() failed to write elements:\n"
                        "Associated writer: \"%s\"\n"
                        "Associated base: \"%s\"\n"
                        "Associated section name: \"%s\"\n%s"),
                      writer->name, base->name, section_name, cg_get_error());
        }

        if (retval == CG_OK)
          retval = cg_elements_partial_write(writer->index,
                                             base->index,
                                             zone_index,
                                             section_index,
                                             elt_start,
                                             elt_end,
                                             _global_connect_s);
        if (retval != CG_OK)
          bft_error(__FILE__, __LINE__, 0,
                    _("cg_elements_partial_write() failed to write elements:\n"
                      "Associated writer: \"%s\"\n"
                      "Associated base: \"%s\"\n"
                      "Associated section name: \"%s\"\n"
                      "Associated range: [%llu, %llu]\n%s\n"),
                    writer->name, base->name, section_name,
                    (unsigned long long) elt_start,
                    (unsigned long long) elt_end,
                    cg_get_error());

        if (sizeof(cgsize_t) > sizeof(cs_gnum_t))
          BFT_FREE(_global_connect_s);
      }

    } /* End if rank == 0 */

    start_id = end_id;

  } /* End of loop on slice advance */

  *global_counter += n_g_sub_elements;

  current_section = current_section->next;

  fvm_gather_slice_destroy(elements_slice);

  BFT_FREE(local_idx);
  BFT_FREE(global_idx_s);
  BFT_FREE(sub_elt_vertex_num);
  if (global_connect_s != global_connect_s_call)
    BFT_FREE(global_connect_s);

  return current_section;
}

#endif /* HAVE_MPI */

/*----------------------------------------------------------------------------
 * Write strided connectivity from tesselated elements to a CGNS file
 * in serial mode.
 *
 * parameters:
 *   export_section   <-- pointer to sections list to export.
 *   writer           <-- pointer to associated writer.
 *   base             <-- pointer to CGNS base structure.
 *   section_id       <-- section identificator number.
 *   global_counter  <-- counter to update the shift after each section export.
 *----------------------------------------------------------------------------*/

static const fvm_writer_section_t *
_export_nodal_tesselated_l(const fvm_writer_section_t  *export_section,
                           const fvm_to_cgns_writer_t  *writer,
                           const fvm_to_cgns_base_t  *base,
                           int  section_id,
                           cs_gnum_t   *global_counter)
{
  int  section_index = -1;
  char  section_name[FVM_CGNS_NAME_SIZE + 1];
  cs_lnum_t   n_sub_elements_max;
  cs_lnum_t   n_buffer_elements_max;
  CGNS_ENUMT(ElementType_t)  cgns_elt_type; /* Definition in cgnslib.h */

  cs_lnum_t   start_id = 0, end_id = 0;
  cs_gnum_t   elt_start = 0, elt_end = 0;

  const cs_lnum_t   *sub_element_idx = NULL;
  cs_lnum_t   *vertex_num = NULL;

  const  int zone_index = 1; /* We always use zone index = 1 */
  const  fvm_writer_section_t *current_section = export_section;
  const  fvm_nodal_section_t *section = current_section->section;
  const  fvm_tesselation_t *tesselation = section->tesselation;
  const  int stride = fvm_nodal_n_vertices_element[current_section->type];

  int  retval = CG_OK;

  assert(current_section != NULL);

  _define_section(current_section, section_id, section_name, &cgns_elt_type);

  sub_element_idx = fvm_tesselation_sub_elt_index(tesselation,
                                                  export_section->type);

  n_buffer_elements_max = section->n_elements;

  fvm_tesselation_get_global_size(section->tesselation,
                                  export_section->type,
                                  NULL,
                                  &n_sub_elements_max);

  if (n_sub_elements_max > n_buffer_elements_max)
    n_buffer_elements_max = n_sub_elements_max;

  BFT_MALLOC(vertex_num, n_buffer_elements_max * stride, cs_lnum_t);

  for (start_id = 0;
       start_id < section->n_elements;
       start_id = end_id) {

    end_id
      = fvm_tesselation_decode(tesselation,
                               current_section->type,
                               start_id,
                               n_buffer_elements_max,
                               export_section->extra_vertex_base,
                               vertex_num);

    /* Print sub-elements connnectivity */

    elt_start = *global_counter + 1 + sub_element_idx[start_id];
    elt_end = *global_counter + sub_element_idx[end_id];

    if (vertex_num != NULL) {
      cgsize_t *_vertex_num = (cgsize_t *)vertex_num;
      if (sizeof(cgsize_t) != sizeof(cs_lnum_t)) {
        int i = 0, n = (elt_end + 1 - elt_start)*stride;
        if (sizeof(cgsize_t) > sizeof(cs_lnum_t))
          BFT_MALLOC(_vertex_num, n, cgsize_t);
        for (i = 0; i < n; i++)
          _vertex_num[i] = vertex_num[i];
      }

      if (start_id == 0) { /* First pass */
        retval = cg_section_partial_write(writer->index,
                                          base->index,
                                          zone_index,
                                          section_name,
                                          cgns_elt_type,
                                          elt_start,
                                          elt_end,
                                          0, /* unsorted boundary elements */
                                          &section_index);
        if (retval != CG_OK)
          bft_error(__FILE__, __LINE__, 0,
                    _("cg_section_partial_write() failed to write "
                      "tesselated elements:\n"
                      "Associated writer: \"%s\"\n"
                      "Associated base: \"%s\"\n"
                      "Associated section name: \"%s\"\n%s"),
                    writer->name, base->name, section_name, cg_get_error());
      }

      if (retval == CG_OK)
        retval = cg_elements_partial_write(writer->index,
                                           base->index,
                                           zone_index,
                                           section_index,
                                           elt_start,
                                           elt_end,
                                           _vertex_num);
      if (retval != CG_OK)
        bft_error(__FILE__, __LINE__, 0,
                  _("cg_elements_partial_write() failed to write elements:\n"
                    "Associated writer: \"%s\"\n"
                    "Associated base: \"%s\"\n"
                    "Associated section name: \"%s\"\n"
                    "Associated range: [%llu, %llu]\n%s\n"),
                  writer->name, base->name, section_name,
                  (unsigned long long) elt_start,
                  (unsigned long long) elt_end,
                  cg_get_error());

      if (sizeof(cgsize_t) > sizeof(cs_lnum_t))
        BFT_FREE(_vertex_num);
    }

  } /* End of loop on parent elements */

  *global_counter += (cs_gnum_t)
    fvm_tesselation_n_sub_elements(tesselation,
                                   current_section->type);

  BFT_FREE(vertex_num);

  return current_section->next;
}

#if defined(HAVE_MPI)
/*----------------------------------------------------------------------------
 * Write polygonal connectivity to a CGNS file in parallel mode
 *
 * parameters:
 *   export_section   <-- pointer to sections list to export.
 *   writer           <-- pointer to associated writer.
 *   mesh             <-- pointer to nodal mesh structure.
 *   base             <-- pointer to CGNS base structure.
 *   section_id       <-- section identificator number.
 *   global_counter   <-- counter to update element shift after each section.
 *   global_s_size    <-- maximum number of entities defined per slice.
 *   global_connect_s_size <-- buffer size to export connectivity
 *   global_connect_s <-- global connectivity array section for elements
 *                        slice global_num_start to global_num_end.
 *                        (output for rank 0, working array only for others)
 *----------------------------------------------------------------------------*/

static const fvm_writer_section_t *
_export_nodal_polygons_g(const fvm_writer_section_t   *export_section,
                         fvm_to_cgns_writer_t   *writer,
                         const fvm_nodal_t      *mesh,
                         fvm_to_cgns_base_t     *base,
                         int                     section_id,
                         cs_gnum_t              *global_counter,
                         cs_gnum_t               global_s_size,
                         cs_gnum_t               global_connect_s_size,
                         cs_gnum_t               global_connect_s_call[])
{
  int  section_index;
  char  section_name[FVM_CGNS_NAME_SIZE + 1];
  cs_lnum_t   elt_id, vertex_id;
  CGNS_ENUMT(ElementType_t) cgns_elt_type; /* Definition in cgnslib.h */

  size_t  connect_end = 0;
  cs_gnum_t   elt_start = 0, elt_end = 0;
  cs_gnum_t   global_num_start = 0, global_num_end = 0;
  cs_lnum_t   *connect_index = NULL;
  cs_lnum_t   *connect_num = NULL;
  cs_gnum_t   *global_connect_s = global_connect_s_call;
  cs_gnum_t   *global_index_s = NULL;
  fvm_gather_slice_t *polygons_slice = NULL;

  const int  zone_index = 1; /* Always = 1 */
  const fvm_writer_section_t *current_section = export_section;
  const fvm_nodal_section_t  *section = current_section->section;
  const cs_lnum_t   n_elements = section->n_elements;
  const cs_lnum_t   *const vertex_index = section->vertex_index;
  const cs_lnum_t   *const vertex_num = section->vertex_num;
  const fvm_io_num_t  *const global_elt_num = section->global_element_num;
  const cs_gnum_t   *vertex_gnum =
    fvm_io_num_get_global_num(mesh->global_vertex_num);

  int  retval = CG_OK;

  assert(current_section != NULL);

  /* Allocate memory for global index array */

  BFT_MALLOC(global_index_s, global_connect_s_size, cs_gnum_t);

  _define_section(current_section, section_id, section_name, &cgns_elt_type);

  /* Create local vertex_index and vertex_num with CGNS standard */

  BFT_MALLOC(connect_index, n_elements + 1, cs_lnum_t);
  BFT_MALLOC(connect_num, vertex_index[n_elements] + n_elements, cs_lnum_t);

  for (elt_id = 0; elt_id < n_elements; elt_id++) {

    connect_index[elt_id] = vertex_index[elt_id] + elt_id;
#if CGNS_VERSION < 3200
    connect_num[connect_end++] = CGNS_ENUMV(NGON_n) + vertex_index[elt_id+1]
                                                    - vertex_index[elt_id];
#else
    connect_num[connect_end++] = vertex_index[elt_id+1] - vertex_index[elt_id];
#endif

    for (vertex_id = vertex_index[elt_id];
         vertex_id < vertex_index[elt_id + 1];
         vertex_id++) {

      if (mesh->global_vertex_num != NULL)
        connect_num[connect_end++] = vertex_gnum[vertex_num[vertex_id] - 1];
      else
        connect_num[connect_end++] = vertex_num[vertex_id];

    }
  }

  connect_index[n_elements] = vertex_index[n_elements] + n_elements;

  polygons_slice
    = fvm_gather_slice_create(global_elt_num,
                              global_s_size,
                              writer->comm);

  while (fvm_gather_slice_advance(polygons_slice,
                                  &global_num_start,
                                  &global_num_end) == 0) {

    /* Gather face->vertices index */

    fvm_gather_slice_index(connect_index,
                           global_index_s,
                           global_elt_num,
                           writer->comm,
                           polygons_slice);

    /* Gather face->vertices connectivity */

    fvm_gather_indexed_numbers(connect_index,
                               connect_num,
                               global_connect_s,
                               NULL, /* already done */
                               global_elt_num,
                               writer->comm,
                               global_index_s,
                               polygons_slice);

    if (writer->rank == 0) {

      elt_start = *global_counter + global_num_start;
      elt_end   = *global_counter + global_num_end - 1;

      if (global_connect_s != NULL) {
        cgsize_t *_global_connect_s = (cgsize_t *)global_connect_s;
        if (sizeof(cgsize_t) != sizeof(cs_gnum_t)) {
          int i = 0, j = 0, n = elt_end + 1 - elt_start;
          for (i = 0; i < n; i++)
#if CGNS_VERSION < 3200
            j += global_connect_s[j] - CGNS_ENUMV(NGON_n) + 1;
#else
            j += global_connect_s[j] + 1;
#endif
          if (sizeof(cgsize_t) > sizeof(cs_gnum_t))
            BFT_MALLOC(_global_connect_s, j, cgsize_t);
          for (i = 0; i < j; i++)
            _global_connect_s[i] = global_connect_s[i];
        }

        if (global_num_start == 1) { /* First pass */
          retval = cg_section_partial_write(writer->index,
                                            base->index,
                                            zone_index,
                                            section_name,
                                            cgns_elt_type,
                                            elt_start,
                                            elt_end,
                                            0, /* unsorted boundary elements */
                                            &section_index);
          if (retval != CG_OK)
            bft_error(__FILE__, __LINE__, 0,
                      _("cg_section_partial_write() failed to write"
                        " polygonal elements:\n"
                        "Associated writer: \"%s\"\n"
                        "Associated base: \"%s\"\n"
                        "Associated section name: \"%s\"\n%s"),
                      writer->name, base->name, section_name, cg_get_error());
        }

        if (retval == CG_OK)
          retval = cg_elements_partial_write(writer->index,
                                             base->index,
                                             zone_index,
                                             section_index,
                                             elt_start,
                                             elt_end,
                                             _global_connect_s);
        if (retval != CG_OK)
          bft_error(__FILE__, __LINE__, 0,
                    _("cg_elements_partial_write() failed to write elements:\n"
                      "Associated writer: \"%s\"\n"
                      "Associated base: \"%s\"\n"
                      "Associated section name: \"%s\"\n"
                      "Associated range: [%llu, %llu]\n%s\n"),
                    writer->name, base->name, section_name,
                    (unsigned long long) elt_start,
                    (unsigned long long) elt_end,
                    cg_get_error());

        if (sizeof(cgsize_t) > sizeof(cs_gnum_t))
          BFT_FREE(_global_connect_s);
      }

    } /* End if rank = 0 */

  } /* End of loop on slice */

  if (elt_end >= elt_start)
    *global_counter += global_num_end - 1;

  /* Free arrays and structure */

  BFT_FREE(connect_index);
  BFT_FREE(connect_num);

  fvm_gather_slice_destroy(polygons_slice);

  BFT_FREE(global_index_s);

  return current_section->next;
}
#endif

/*----------------------------------------------------------------------------
 * Write polygonal connectivity to a CGNS file in serial mode
 *
 * parameters:
 *   export_section <-- pointer to sections list to export.
 *   writer         <-- pointer to associated writer.
 *   base           <-- pointer to CGNS base structure.
 *   section_id       <-- section identificator number.
 *   global_counter <-- counter to update the shift after each section export.
 *----------------------------------------------------------------------------*/

static const fvm_writer_section_t *
_export_nodal_polygons_l(const fvm_writer_section_t  *export_section,
                         const fvm_to_cgns_writer_t  *writer,
                         const fvm_to_cgns_base_t    *base,
                         int                          section_id,
                         cs_gnum_t                   *global_counter)
{
  int  i, j, section_index;
  char  section_name[FVM_CGNS_NAME_SIZE + 1];
  CGNS_ENUMT(ElementType_t)  cgns_elt_type; /* Definition in cgnslib.h */

  cs_lnum_t   connect_size = 0;
  cs_gnum_t   elt_start = 0, elt_end = 0;
  cgsize_t  *connect = NULL;

  const  int  zone_index = 1; /* We always use zone index = 1 */
  const fvm_writer_section_t *current_section = export_section;
  const fvm_nodal_section_t *section = current_section->section;

  int  retval = CG_OK;

  assert(current_section != NULL);

  _define_section(export_section, section_id, section_name, &cgns_elt_type);

  elt_start = *global_counter + 1;
  elt_end = *global_counter + section->n_elements;

  BFT_MALLOC(connect,
             section->n_elements + section->connectivity_size,
             cgsize_t);

  for (j = 0; j < section->n_elements; j++) {
#if CGNS_VERSION < 3200
    connect[connect_size++] = CGNS_ENUMV(NGON_n) + section->vertex_index[j+1]
                                                 - section->vertex_index[j];
#else
    connect[connect_size++] =    section->vertex_index[j+1]
                               - section->vertex_index[j];
#endif

    for (i = section->vertex_index[j]; i < section->vertex_index[j+1]; i++)
      connect[connect_size++] = section->vertex_num[i];
  }

  if (connect != NULL)
    retval = cg_section_write(writer->index,
                              base->index,
                              zone_index,
                              section_name,
                              cgns_elt_type,
                              elt_start,
                              elt_end,
                              0, /* unsorted boundary elements */
                              connect,
                              &section_index);

  if (retval != CG_OK)
    bft_error(__FILE__, __LINE__, 0,
              _("cg_section_write() failed to write polygonal elements:\n"
                "Associated writer: \"%s\"\n"
                "Associated base: \"%s\"\n"
                "Associated section name: \"%s\"\n%s"),
              writer->name, base->name, section_name, cg_get_error());

  *global_counter += section->n_elements;

  BFT_FREE(connect);

  return current_section->next;
}

/*----------------------------------------------------------------------------
 * Write a per element solution field to a CGNS file
 *
 * parameters:
 *   export_list        <-- pointer to section helper structure
 *   helper             <-- pointer to general writer helper structure
 *   writer             <-- pointer to associated writer.
 *   base               <-- pointer to CGNS base structure.
 *   fieldlabel         <-- variable name.
 *   solution_index     <-- index of the associated CGNS solution structure.
 *   input_dim          <-- input field dimension.
 *   cgns_datatype      <-- indicates the CGNS data type of output field.
 *   interlace          <-- indicates if variable in memory is interlaced.
 *   n_parent_lists     <-- indicates if variable values are to be obtained
 *                          directly through the local entity index (when 0) or
 *                          through the parent entity numbers (when 1 or more).
 *   parent_num_shift   <-- parent number to value array index shifts
 *                          size: n_parent_lists.
 *   datatype           <-- indicates the data type of (source) field values.
 *   field_values       <-- array of associated field value arrays.
 *   output_buffer_size <-- output buffer size
 *   output_buffer      --- output buffer (size output_buffer_size)
 *----------------------------------------------------------------------------*/

static void
_export_field_e(const fvm_writer_section_t      *export_list,
                fvm_writer_field_helper_t       *helper,
                const fvm_to_cgns_writer_t      *writer,
                const fvm_to_cgns_base_t        *base,
                const char                      *fieldlabel,
                int                              solution_index,
                int                              input_dim,
                CGNS_ENUMT(DataType_t)           cgns_datatype,
                cs_interlace_t                   interlace,
                int                              n_parent_lists,
                const cs_lnum_t                  parent_num_shift[],
                cs_datatype_t                    datatype,
                const void                *const field_values[],
                size_t                           output_buffer_size,
                void                      *const output_buffer)
{
  int  i;
  size_t  output_size;

  const fvm_writer_section_t  *current_section = NULL;

  int output_dim = fvm_writer_field_helper_field_dim(helper);

  /* Loop on dimension (always de-interlace vectors) */

  for (i = 0 ; i < output_dim ; i++) {

    cgsize_t partial_write_idx_start = 1;

    for (current_section = export_list;
         current_section != NULL;
         current_section = current_section->next) {

      while (fvm_writer_field_helper_step_e(helper,
                                            current_section,
                                            input_dim,
                                            i,
                                            interlace,
                                            n_parent_lists,
                                            parent_num_shift,
                                            datatype,
                                            field_values,
                                            output_buffer,
                                            output_buffer_size,
                                            &output_size) == 0) {

        if (writer->rank == 0) {

          int field_index;
          int retval = CG_OK;
          cgsize_t partial_write_idx_end
            = partial_write_idx_start + output_size - 1;

          const int  shift = FVM_CGNS_NAME_SIZE + 1;
          const int  zone_index = 1; /* We always work on index zone = 1 */

          retval = cg_field_partial_write(writer->index,
                                          base->index,
                                          zone_index,
                                          solution_index,
                                          cgns_datatype,
                                          &fieldlabel[i * shift],
                                          &partial_write_idx_start,
                                          &partial_write_idx_end,
                                          output_buffer,
                                          &field_index);

          if (retval != CG_OK)
            bft_error(__FILE__, __LINE__, 0,
                      _("cg_field_partial_write() failed to write "
                        "field values:\n\"%s\"\n"
                        "Associated writer: \"%s\"\n"
                        "Associated base: \"%s\"\n%s"),
                      fieldlabel + i*shift, writer->name, base->name,
                      cg_get_error());

          partial_write_idx_start = partial_write_idx_end + 1;

        } /* end of write for rank 0 */

      } /* end of field writer helper step */

    } /* end of loop on sections */

  } /* end of loop on spatial dimension */

}

/*----------------------------------------------------------------------------
 * Write a per node solution field to a CGNS file
 *
 * parameters:
 *   mesh               <-- pointer to nodal mesh  structure that should
 *                          be written.
 *   helper             <-- pointer to general writer helper structure
 *   writer             <-- pointer to associated writer.
 *   base               <-- pointer to CGNS base structure.
 *   fieldlabel         <-- variable name.
 *   solution_index     <-- index of the associated CGNS solution structure.
 *   input_dim          <-- input field dimension.
 *   cgns_datatype      <-- indicates the CGNS data type of output field.
 *   interlace          <-- indicates if variable in memory is interlaced.
 *   n_parent_lists     <-- indicates if variable values are to be obtained
 *                          directly through the local entity index (when 0) or
 *                          through the parent entity numbers (when 1 or more).
 *   parent_num_shift   <-- parent number to value array index shifts
 *                          size: n_parent_lists.
 *   datatype           <-- indicates the data type of (source) field values.
 *   field_values       <-- array of associated field value arrays.
 *   output_buffer_size <-- output buffer size
 *   output_buffer      --- output buffer (size output_buffer_size)
 *----------------------------------------------------------------------------*/

static void
_export_field_n(const fvm_nodal_t               *mesh,
                fvm_writer_field_helper_t       *helper,
                const fvm_to_cgns_writer_t      *writer,
                const fvm_to_cgns_base_t        *base,
                const char                      *fieldlabel,
                int                              solution_index,
                int                              input_dim,
                CGNS_ENUMT(DataType_t)           cgns_datatype,
                cs_interlace_t                   interlace,
                int                              n_parent_lists,
                const cs_lnum_t                  parent_num_shift[],
                cs_datatype_t                    datatype,
                const void                *const field_values[],
                size_t                           output_buffer_size,
                void                      *const output_buffer)
{
  int  i;
  size_t  output_size;

  int output_dim = fvm_writer_field_helper_field_dim(helper);

  /* Loop on dimension (always de-interlace vectors) */

  for (i = 0 ; i < output_dim ; i++) {

    cgsize_t partial_write_idx_start = 1;

    while (fvm_writer_field_helper_step_n(helper,
                                          mesh,
                                          input_dim,
                                          i,
                                          interlace,
                                          n_parent_lists,
                                          parent_num_shift,
                                          datatype,
                                          field_values,
                                          output_buffer,
                                          output_buffer_size,
                                          &output_size) == 0) {

      if (writer->rank == 0) {

        int field_index;
        int retval = CG_OK;
        cgsize_t partial_write_idx_end
          = partial_write_idx_start + output_size - 1;

        const int  shift = FVM_CGNS_NAME_SIZE + 1;
        const int  zone_index = 1; /* We always work on index zone = 1 */

        retval = cg_field_partial_write(writer->index,
                                        base->index,
                                        zone_index,
                                        solution_index,
                                        cgns_datatype,
                                        &fieldlabel[i * shift],
                                        &partial_write_idx_start,
                                        &partial_write_idx_end,
                                        output_buffer,
                                        &field_index);

        if (retval != CG_OK)
          bft_error(__FILE__, __LINE__, 0,
                    _("cg_field_partial_write() failed to write "
                      "field values:\n\"%s\"\n"
                      "Associated writer: \"%s\"\n"
                      "Associated base: \"%s\"\n%s"),
                    fieldlabel + i*shift, writer->name, base->name,
                    cg_get_error());

        partial_write_idx_start = partial_write_idx_end + 1;

      } /* end of write for rank 0 */

    } /* end of field writer helper step */

  } /* end of loop on spatial dimension */

}

/*----------------------------------------------------------------------------
 * Create time-dependent data structure in CGNS file: Base Iterative Data
 * structure and Zone Iterative Data structure
 *
 * parameters:
 *   writer  <-- CGNS writer structure
 *----------------------------------------------------------------------------*/

static void
_create_timedependent_data(fvm_to_cgns_writer_t  *writer)
{
  int     base_id, j, name_len;
  cgsize_t dim[2];

  double *time_values = NULL;
  int    *time_steps = NULL;
  char   *sol_names = NULL;

  const int  zone_index = 1;

  int     retval = CG_OK;
  int     sol_id = -1;

  assert(writer->bases != NULL);
  assert(writer->is_open == true);

  /* Create structures for time-dependent data */
  /*-------------------------------------------*/

  /* Create a BaseIterativeData */

  for (base_id = 0; base_id < writer->n_bases; base_id++) {
    fvm_to_cgns_base_t *base = writer->bases[base_id];

    if (base->n_sols == 0)
      continue;

    retval = cg_biter_write(writer->index,
                            base->index,
                            "BaseIterativeData_t",
                            base->n_sols);

    if (retval != CG_OK)
      bft_error(__FILE__, __LINE__, 0,
                _("cg_biter_write() failed to create a BaseIterativeData\n"
                  "Associated writer:\"%s\" :\n"
                  "Associated base:\"%s\"\n%s"),
                writer->filename, base->name,cg_get_error());

    retval = cg_goto(writer->index,
                     base->index,
                     "BaseIterativeData_t",
                     1,
                     "end");

    if (retval == CG_OK) {

      BFT_MALLOC(time_values, base->n_sols, double);
      BFT_MALLOC(time_steps, base->n_sols, int);

      for (sol_id = 0; sol_id < base->n_sols; sol_id++) {
        time_values[sol_id] = base->solutions[sol_id]->time_value;
        time_steps[sol_id] = base->solutions[sol_id]->time_step;
      }

      dim[0] = sol_id;
      retval = cg_array_write("TimeValues",
                              CGNS_ENUMV(RealDouble),
                              1,
                              dim,
                              time_values);

      if (retval != CG_OK)
        bft_error(__FILE__, __LINE__, 0,
                  _("cg_array_write() failed to write TimeValues\n"
                    "Associated writer:\"%s\" :\n"
                    "Associated base:\"%s\"\n%s"),
                  writer->filename, base->name,cg_get_error());

      dim[0] = sol_id;
      retval = cg_array_write("IterationValues",
                              CGNS_ENUMV(Integer),
                              1,
                              dim,
                              time_steps);

      if (retval != CG_OK)
        bft_error(__FILE__, __LINE__, 0,
                  _("cg_array_write failed to write IterationValues\n"
                    "Associated writer:\"%s\" :\n"
                    "Associated base:\"%s\"\n%s"),
                  writer->filename, base->name,cg_get_error());

      BFT_FREE(time_values);
      BFT_FREE(time_steps);
    }

    /* Create a ZoneIterativeData */

    retval = cg_ziter_write(writer->index,
                            base->index,
                            zone_index,
                            "ZoneIterativeData");

    if (retval != CG_OK)
      bft_error(__FILE__, __LINE__, 0,
                _("cg_ziter_write() failed to create a ZoneIterativeData\n"
                  "Associated writer:\"%s\" :\n"
                  "Associated base:\"%s\"\n%s"),
                writer->filename, base->name,cg_get_error());

    retval = cg_goto(writer->index,
                     base->index,
                     "Zone_t", zone_index,
                     "ZoneIterativeData_t", 1,
                     "end");

    if (retval == CG_OK) {

      /* Write solution names in a array */

      dim[0] = FVM_CGNS_NAME_SIZE;
      dim[1] = sol_id;

      BFT_MALLOC(sol_names, dim[0] * dim[1] , char);

      for (j = 0; j < dim[0] * dim[1]; j++)
        sol_names[j] = ' ';

      for (sol_id = 0; sol_id < base->n_sols; sol_id++) {
        name_len = strlen(base->solutions[sol_id]->name);
        strncpy(sol_names + sol_id * FVM_CGNS_NAME_SIZE,
                base->solutions[sol_id]->name, name_len);
      }

      retval = cg_array_write("FlowSolutionPointers",
                              CGNS_ENUMV(Character),
                              2,
                              dim,
                              sol_names);

      if (retval != CG_OK)
        bft_error(__FILE__, __LINE__, 0,
                  _("cg_array_write() failed to write FlowSolutionPointers\n"
                    "Associated writer:\"%s\" :\n"
                    "Associated base:\"%s\"\n%s"),
                  writer->filename, base->name, cg_get_error());

      BFT_FREE(sol_names);

    }

    retval = cg_simulation_type_write(writer->index,
                                      base->index,
                                      CGNS_ENUMV(TimeAccurate));

    if (retval != CG_OK)
      bft_error(__FILE__, __LINE__, 0,
                _("cg_simulation_type_write() failed\n"
                  "Associated writer:\"%s\" :\n"
                  "Associated base:\"%s\"\n%s"),
                writer->filename, base->name, cg_get_error());

  } /* End of loop on bases */

}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Returns number of library version strings associated with the CGNS format.
 *
 * returns:
 *   number of library version strings associated with the CGNS format.
 *----------------------------------------------------------------------------*/

int
fvm_to_cgns_n_version_strings(void)
{
  return 1;
}

/*----------------------------------------------------------------------------
 * Returns a library version string associated with the CGNS format.
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
fvm_to_cgns_version_string(int string_index,
                           int compile_time_version)
{
  const char * retval = NULL;

  if (string_index == 0) {
    snprintf(_cgns_version_string, 31, "CGNS %d.%d.%d\n",
             CGNS_VERSION/1000,
             (CGNS_VERSION % 1000) / 100,
             (CGNS_VERSION % 100) / 10);
    _cgns_version_string[31] = '\0';
    retval = _cgns_version_string;
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Initialize FVM to CGNS file writer.
 *
 * Options are:
 *   discard_polygons    do not output polygons or related values
 *   discard_polyhedra   do not output polyhedra or related values
 *   divide_polygons     tesselate polygons with triangles
 *   adf                 use ADF file type
 *   hdf5                use HDF5 file type (default if available)
 *
 * As CGNS does not handle polyhedral elements, polyhedra are automatically
 * tesselated with tetrahedra and pyramids (adding a vertex near each
 * polyhedron's center) unless discarded.
 *
 * parameters:
 *   name           <-- base output case name.
 *   options        <-- whitespace separated, lowercase options list
 *   time_dependecy <-- indicates if and how meshes will change with time
 *   comm           <-- associated MPI communicator.
 *
 * returns:
 *   pointer to opaque CGNS writer structure.
 *----------------------------------------------------------------------------*/

#if defined(HAVE_MPI)
void *
fvm_to_cgns_init_writer(const char             *name,
                        const char             *path,
                        const char             *options,
                        fvm_writer_time_dep_t   time_dependency,
                        MPI_Comm                comm)
#else
void *
fvm_to_cgns_init_writer(const char             *name,
                        const char             *path,
                        const char             *options,
                        fvm_writer_time_dep_t   time_dependency)
#endif
{
  int  i, writer_index;
  int  filename_length, name_length, path_length;
  bool force_adf = false, force_hdf5 = false;

  fvm_to_cgns_writer_t  *writer = NULL;

  /* Initialize writer */

  BFT_MALLOC(writer, 1, fvm_to_cgns_writer_t);

  /* Mesh time dependency */

  writer->time_dependency = time_dependency;

  /* Writer name */

  name_length = strlen(name);
  if (name_length == 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Empty CGNS filename."));
  BFT_MALLOC(writer->name, name_length + 1, char);
  strcpy(writer->name, name);

  for (i = 0; i < name_length; i++) {
    if (writer->name[i] == ' ' || writer->name[i] == '\t')
      writer->name[i] = '_';
  }

  /* Writer's associated filename */

  if (path != NULL)
    path_length = strlen(path);
  else
    path_length = 0;
  filename_length = path_length + name_length + strlen(".cgns") + 1;
  BFT_MALLOC(writer->filename, filename_length, char);

  if (path != NULL)
    strcpy(writer->filename, path);
  else
    writer->filename[0] = '\0';

  strcat(writer->filename, writer->name);
  strcat(writer->filename, ".cgns");

  /* CGNS Base structure */

  writer->n_bases = 0;
  writer->bases = NULL;

  /* Other variables */

  writer->n_time_steps = 0;
  writer->time_steps = NULL;
  writer->time_values = NULL;
  writer->rank = 0;
  writer->n_ranks = 1;

  writer->discard_polygons = false;
  writer->discard_polyhedra = false;
  writer->divide_polygons = false;

/* As CGNS does not handle polyhedral elements, polyhedra are automatically
 * tesselated with tetrahedra and pyramids (adding a vertex near each
 * polyhedron's center) unless discarded. */

#if defined(HAVE_MPI)
  {
    int mpi_flag, rank, n_ranks;
    MPI_Initialized(&mpi_flag);

    if (mpi_flag && comm != MPI_COMM_NULL) {
      writer->comm = comm;
      MPI_Comm_rank(writer->comm, &rank);
      MPI_Comm_size(writer->comm, &n_ranks);
      writer->rank = rank;
      writer->n_ranks = n_ranks;
    }
    else
      writer->comm = MPI_COMM_NULL;
  }
#endif /* defined(HAVE_MPI) */

  /* Parse options */

  if (options != NULL) {
    int i1, i2, l_opt;
    int l_tot = strlen(options);

    i1 = 0; i2 = 0;
    while (i1 < l_tot) {
      for (i2 = i1; i2 < l_tot && options[i2] != ' '; i2++);

      l_opt = i2 - i1;

      if (   (l_opt == 16)
          && (strncmp(options + i1, "discard_polygons", l_opt) == 0))
        writer->discard_polygons = true;

      else if (   (l_opt == 17)
               && (strncmp(options + i1, "discard_polyhedra", l_opt) == 0))
        writer->discard_polyhedra = true;

      else if (   (l_opt == 15)
               && (strncmp(options + i1, "divide_polygons", l_opt) == 0))
        writer->divide_polygons = true;

      else if (   (l_opt == 3)
               && (strncmp(options + i1, "adf", l_opt) == 0))
        force_adf = false;

      else if (   (l_opt == 4)
               && (strncmp(options + i1, "hdf5", l_opt) == 0))
        force_hdf5 = false;

      for (i1 = i2 + 1; i1 < l_tot && options[i1] == ' '; i1++);
    }
  }

  /* Open CNGS file */

  writer->is_open = false;

  if (writer->rank == 0) {

    if (force_adf == true)
      cg_set_file_type(CG_FILE_ADF);

    if (force_hdf5 == true)
      cg_set_file_type(CG_FILE_HDF5);

    if (cg_open(writer->filename, CG_MODE_WRITE, &writer_index) != CG_OK)
      bft_error(__FILE__, __LINE__, 0,
                _("cg_open() failed to open file \"%s\" : \n%s"),
                writer->filename, cg_get_error());

    writer->is_open = true;

  }

#if defined(HAVE_MPI)
  if (writer->n_ranks > 1)
    MPI_Bcast(&writer_index, 1, MPI_INT, 0, writer->comm);
#endif

  writer->index = writer_index;

  return writer;
}

/*----------------------------------------------------------------------------
 * Finalize FVM to CGNS file writer.
 *
 * parameters:
 *   this_writer_p <-- pointer to opaque CGNS writer structure.
 *
 * returns:
 *   NULL pointer.
 *----------------------------------------------------------------------------*/

void *
fvm_to_cgns_finalize_writer(void  *this_writer_p)
{
  int i;

  fvm_to_cgns_writer_t  *writer
                        = (fvm_to_cgns_writer_t *)this_writer_p;

  assert(writer != NULL);

  if (writer->rank == 0) {

    if (writer->bases != NULL) {

      /* Create index for time-dependent data */
      /*--------------------------------------*/

      _create_timedependent_data(writer);

    }

    /* Close CGNS File */
    /*-----------------*/

    if (writer->is_open == true) {

      if (cg_close(writer->index) != CG_OK)
        bft_error(__FILE__, __LINE__, 0,
                  _("cg_close() failed to close file \"%s\" :\n%s"),
                  writer->filename, cg_get_error());

    }

  } /* End if rank = 0 */

  /* Free structures */
  /*-----------------*/

  /* Free names */

  BFT_FREE(writer->name);
  BFT_FREE(writer->filename);
  BFT_FREE(writer->time_values);
  BFT_FREE(writer->time_steps);

  /* Free fvm_to_cgns_base structure */

  for (i = 0; i < writer->n_bases; i++)
    writer->bases[i] = _del_base(writer->bases[i]);

  BFT_FREE(writer->bases);

  /* Free fvm_to_cgns_writer structure */

  BFT_FREE(writer);
  return NULL;
}

/*----------------------------------------------------------------------------
 * Associate new time value with a writer structure if necessary.
 *
 * parameters:
 *   this_writer_p <-- pointer to associated writer
 *   time_step     <-- time step number
 *   time_value    <-- time_value number
 *----------------------------------------------------------------------------*/

void
fvm_to_cgns_set_mesh_time(void     *this_writer_p,
                          int       time_step,
                          double    time_value)
{
  int n_vals;

  fvm_to_cgns_writer_t  *writer
                        = (fvm_to_cgns_writer_t *)this_writer_p;

  static char time_value_err_string[] =
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

}

/*----------------------------------------------------------------------------
 * Indicate if elements of a given type in a mesh associated with a given
 * CGNS file writer need to be tesselated.
 *
 * parameters:
 *   this_writer_p <-- pointer to associated writer
 *   mesh          <-- pointer to nodal mesh structure that should be written
 *   element_type  <-- element type we are interested in
 *
 * returns:
 *   1 if tesselation of the given element type is needed, 0 otherwise
 *----------------------------------------------------------------------------*/

int
fvm_to_cgns_needs_tesselation(fvm_writer_t       *this_writer_p,
                              const fvm_nodal_t  *mesh,
                              fvm_element_t       element_type)
{
  int  i;
  int  retval = 0;
  fvm_to_cgns_writer_t  *this_writer
                             = (fvm_to_cgns_writer_t *)this_writer_p;

  const int  export_dim = fvm_nodal_get_max_entity_dim(mesh);

  /* Initial count and allocation */

  if (   (   element_type == FVM_FACE_POLY
          && this_writer->divide_polygons == true)
      || (   element_type == FVM_CELL_POLY
          && this_writer->discard_polyhedra == false)) {

    for (i = 0 ; i < mesh->n_sections ; i++) {

      const fvm_nodal_section_t  *const  section = mesh->sections[i];

      /* Output if entity dimension equal to highest in mesh
         (i.e. no output of faces if cells present, or edges
         if cells or faces) */

      if (section->entity_dim == export_dim) {
        if (section->type == element_type)
          retval = 1;
      }

    }

  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Write nodal mesh to a CGNS file
 *
 * parameters:
 *   this_writer_p <-- pointer to associated writer.
 *   mesh          <-- pointer to nodal mesh structure that should be written.
 *----------------------------------------------------------------------------*/

void
fvm_to_cgns_export_nodal(void               *this_writer_p,
                         const fvm_nodal_t  *mesh)
{
  char  base_name[FVM_CGNS_NAME_SIZE+1];
  int   base_index;

  int section_id = 0;
  cs_gnum_t   n_g_vertices = 0;
  cs_gnum_t   global_counter = 0;
  cs_gnum_t   global_connect_s_size = 0;
  cs_gnum_t   global_s_size = 0;

  const fvm_writer_section_t  *export_section = NULL;
  fvm_writer_section_t  *export_list = NULL;
  fvm_to_cgns_base_t *base = NULL;
  fvm_to_cgns_writer_t  *writer
                        = (fvm_to_cgns_writer_t *)this_writer_p;

  cs_gnum_t   *global_connect_s = NULL;

  const int   n_ranks = writer->n_ranks;


  /* Initialization */
  /*----------------*/

  /* Clean mesh->name */

  strncpy(base_name, mesh->name, FVM_CGNS_NAME_SIZE);
  base_name[FVM_CGNS_NAME_SIZE] = '\0';

  /* Get CGNS base index */

  base_index = _get_base_index(writer,
                               base_name);

  if (base_index == 0 )
    base_index = _add_base(writer,
                           base_name,
                           mesh);

  base = writer->bases[base_index - 1];

  /* Buffer and global sizes required in parallel mode */

  fvm_writer_def_nodal_buf_size(mesh,
                                n_ranks,
                                10,
                                5,
                                &n_g_vertices,
                                NULL,
                                &global_s_size,
                                &global_connect_s_size);

  /* Avoid too many small communications with large processor counts */

  if (n_ranks > 1 && global_connect_s_size > 0) {
    size_t min_buffer_size =   cs_parall_get_min_coll_buf_size()
                             / sizeof(cs_gnum_t);
    if (min_buffer_size > global_connect_s_size) {
      cs_gnum_t global_s_size_min = global_s_size * (  min_buffer_size
                                                      / global_connect_s_size);
      if (global_s_size_min > global_s_size) {
        global_s_size = global_s_size_min;
        global_connect_s_size = min_buffer_size;
      }
    }
  }

  /* Build list of sections that are used here, in order of output */

  export_list = fvm_writer_export_list(mesh,
                                       fvm_nodal_get_max_entity_dim(mesh),
                                       true,
                                       writer->discard_polygons,
                                       writer->discard_polyhedra,
                                       writer->divide_polygons,
                                       true);

  /* Create a zone */
  /*---------------*/

  _add_zone(mesh,
            writer,
            base,
            export_list,
            n_g_vertices);

  /* Vertex coordinates */
  /*--------------------*/

#if defined(HAVE_MPI)

  if (n_ranks > 1)
    _export_vertex_coords_g(writer,
                            mesh,
                            base,
                            global_s_size);

#endif /* HAVE_MPI */

  if (n_ranks == 1)
    _export_vertex_coords_l(writer,
                            mesh,
                            base);

  /* Element connectivity */
  /*----------------------*/

  if (n_ranks > 1)
    BFT_MALLOC(global_connect_s, global_connect_s_size, cs_gnum_t);

  export_section = export_list;

  while (export_section != NULL) {

    const fvm_nodal_section_t  *section = export_section->section;

    /* update section_id (used in section name) */

    section_id++;

    /* Output for strided (regular) element types */
    /*--------------------------------------------*/

    if (section->stride > 0) {

#if defined(HAVE_MPI)
      if (n_ranks > 1)
        export_section = _export_connect_g(export_section,
                                           writer,
                                           mesh,
                                           base,
                                           section_id,
                                           &global_counter,
                                           global_s_size,
                                           global_connect_s);
#endif

      if (n_ranks == 1)
        export_section = _export_connect_l(export_section,
                                           writer,
                                           base,
                                           section_id,
                                           &global_counter);

    }

    /* Output for tesselated polygons or polyhedra */
    /*---------------------------------------------*/

    else if (export_section->type != section->type) {

#if defined(HAVE_MPI)
      if (n_ranks > 1)
        export_section = _export_nodal_tesselated_g(export_section,
                                                    writer,
                                                    mesh,
                                                    base,
                                                    section_id,
                                                    &global_counter,
                                                    global_s_size,
                                                    global_connect_s_size,
                                                    global_connect_s);
#endif

      if (n_ranks == 1)
        export_section = _export_nodal_tesselated_l(export_section,
                                                    writer,
                                                    base,
                                                    section_id,
                                                    &global_counter);

    }

    /* Output for polygons */
    /*---------------------*/

    else if (export_section->type == FVM_FACE_POLY) {

#if defined(HAVE_MPI)
      if (n_ranks > 1)
        export_section = _export_nodal_polygons_g(export_section,
                                                  writer,
                                                  mesh,
                                                  base,
                                                  section_id,
                                                  &global_counter,
                                                  global_s_size,
                                                  global_connect_s_size,
                                                  global_connect_s);
#endif

      if (n_ranks == 1)
        export_section =  _export_nodal_polygons_l(export_section,
                                                   writer,
                                                   base,
                                                   section_id,
                                                   &global_counter);

    }

  } /* End of loop on sections */

  /* Free buffers */

  BFT_FREE(export_list);

#if defined(HAVE_MPI)
  if (n_ranks > 1) {
    BFT_FREE(global_connect_s);
  }
#endif

}

/*----------------------------------------------------------------------------
 * Write field associated with a nodal mesh to a CGNS file.
 *
 * Assigning a negative value to the time step indicates a time-independent
 * field (in which case the time_value argument is unused).
 *
 * parameters:
 *   this_writer_p    <-- pointer to associated writer
 *   mesh             <-- pointer to associated nodal mesh structure
 *   name             <-- variable name
 *   location         <-- fvm grid location (nodes or elements)
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
fvm_to_cgns_export_field(void                   *this_writer_p,
                         const fvm_nodal_t      *mesh,
                         const char             *name,
                         fvm_writer_var_loc_t    location,
                         int                     dimension,
                         cs_interlace_t          interlace,
                         int                     n_parent_lists,
                         const cs_lnum_t         parent_num_shift[],
                         cs_datatype_t           datatype,
                         int                     time_step,
                         double                  time_value,
                         const void       *const field_values[])
{
  char   base_name[FVM_CGNS_NAME_SIZE+1];
  char   field_name[FVM_CGNS_NAME_SIZE+1];
  cs_datatype_t  export_datatype = CS_DATATYPE_NULL;
  int    output_dim;

  size_t  input_size = 0;
  size_t  output_size = 0;
  size_t  min_var_buffer_size = 0;
  size_t  var_buffer_size = 0;
  double *var_buffer = NULL;

  fvm_writer_field_helper_t  *helper = NULL;
  fvm_writer_section_t  *export_list = NULL;

  fvm_to_cgns_solution_t *solution;

  char  *field_label = NULL;
  fvm_to_cgns_writer_t  *writer
                        = (fvm_to_cgns_writer_t *)this_writer_p;
  int base_index = 0;
  int sol_index = 0;
  CGNS_ENUMT(DataType_t)      cgns_datatype = CGNS_ENUMV(DataTypeNull);
  CGNS_ENUMT(GridLocation_t)  cgns_location = CGNS_ENUMV(GridLocationNull);

  const int  rank = writer->rank;
  const int  n_ranks = writer->n_ranks;

  /* Initialization */
  /*----------------*/

  /* FVM datatype conversion to CGNS */

  switch (datatype) {
  case CS_DOUBLE:
    cgns_datatype = CGNS_ENUMV(RealDouble);
    export_datatype = CS_DOUBLE;
    break;
  case CS_FLOAT:
    cgns_datatype = CGNS_ENUMV(RealSingle);
    export_datatype = CS_FLOAT;
    break;
  case CS_UINT32:
    cgns_datatype = CGNS_ENUMV(Integer);
    export_datatype = CS_INT32;
    break;
  case CS_UINT64:
    cgns_datatype = CGNS_ENUMV(Integer);
    export_datatype = CS_INT32;
    break;
  case CS_INT32:
    cgns_datatype = CGNS_ENUMV(Integer);
    export_datatype = CS_INT32;
    break;
  case CS_INT64:
    cgns_datatype = CGNS_ENUMV(Integer);
    export_datatype = CS_INT32;
    break;
  default:
    assert(0);
  }

  /* Set CGNS location */

  if (location == FVM_WRITER_PER_NODE)
    cgns_location = CGNS_ENUMV(Vertex);
  else if (location == FVM_WRITER_PER_ELEMENT)
    cgns_location = CGNS_ENUMV(CellCenter);

  /* Cell dimension */

  output_dim = dimension;

  assert(output_dim > 0);

  if (dimension == 2)
    output_dim = 3;
  else if (dimension > 3 && dimension != 6 && dimension != 9)
    bft_error(__FILE__, __LINE__, 0,
              _("Data of dimension %d not handled"), dimension);

  /* Get CGNS base index */
  /*---------------------*/

  strncpy(base_name, mesh->name, FVM_CGNS_NAME_SIZE);
  base_name[FVM_CGNS_NAME_SIZE] = '\0';

  base_index = _get_base_index(writer, base_name);

  if (base_index == 0 )
    base_index = _add_base(writer,
                           base_name,
                           mesh);

  /* Get CGNS solution index */
  /*-------------------------*/

  sol_index = _get_solution_index(writer,
                                  base_index,
                                  time_step,
                                  time_value,
                                  cgns_location);

  if (sol_index == 0)
    sol_index = _add_solution(writer,
                              base_index,
                              time_step,
                              time_value,
                              cgns_location);

  solution = writer->bases[base_index - 1]->solutions[sol_index - 1];
  assert(solution->location == cgns_location);

  /* Field_Name adaptation if necessary */
  /*-----------------------------------*/

  if (rank == 0) {

    int        i, shift, pos;
    char      *tmp;

    strncpy(field_name, name, FVM_CGNS_NAME_SIZE);
    field_name[FVM_CGNS_NAME_SIZE] = '\0';

    shift = FVM_CGNS_NAME_SIZE + 1;
    BFT_MALLOC(field_label, output_dim * shift, char);

    for (pos = strlen(field_name) - 1;
         pos > 0 && (field_name[pos] == ' ' || field_name[pos] == '\t');
         pos--);
    pos++;

    for (i = 0; i < output_dim; i++) {

      tmp = field_label + (i * shift);
      strncpy(tmp, field_name, shift - 1);
      tmp[shift - 1] = '\0';

      if (output_dim > 1) {

        if (output_dim == 3) {

          const char *comp[] = {"X", "Y", "Z"};

          if (pos > shift - 2)
            pos = shift - 2;

          tmp[pos    ] = comp[i][0];
          tmp[pos + 1] = '\0';

        }
        else if (output_dim == 6) {

          const char *comp[] = {"XX", "YY", "ZZ", "XY", "XZ", "YZ"};

          if (pos > shift - 3)
            pos = shift - 3;

          tmp[pos    ] = comp[i][0];
          tmp[pos + 1] = comp[i][1];
          tmp[pos + 2] = '\0';

        }
        else if (output_dim == 9) {

          const char *comp[] = {"XX", "XY", "XZ",
                                "YX", "YY", "YZ",
                                "ZX", "ZY", "ZZ"};

          if (pos > shift - 3)
            pos = shift - 3;

          tmp[pos    ] = comp[i][0];
          tmp[pos + 1] = comp[i][1];
          tmp[pos + 2] = '\0';

        }
      }
    } /* End of loop on output dimension */

  } /* End if rank == 0 */

  /* Initialize writer helper */
  /*--------------------------*/

  export_list = fvm_writer_export_list(mesh,
                                       fvm_nodal_get_max_entity_dim(mesh),
                                       true,
                                       writer->discard_polygons,
                                       writer->discard_polyhedra,
                                       writer->divide_polygons,
                                       true);

  helper = fvm_writer_field_helper_create(mesh,
                                          export_list,
                                          output_dim,
                                          CS_NO_INTERLACE,
                                          export_datatype,
                                          location);

#if defined(HAVE_MPI)

  fvm_writer_field_helper_init_g(helper,
                                 export_list,
                                 mesh,
                                 writer->comm);

#endif

  /* Buffer size computation and allocation */

  fvm_writer_field_helper_get_size(helper,
                                   &input_size,
                                   &output_size,
                                   NULL,
                                   &min_var_buffer_size);

  /* Slicing allows for arbitrary buffer size, but should be small enough
     to add little additional memory requirement (in proportion), large
     enough to limit number of write and gather calls. */

  if (n_ranks > 1)
    var_buffer_size = input_size / n_ranks;
  else
    var_buffer_size = input_size / 4;

  var_buffer_size = CS_MAX(var_buffer_size, min_var_buffer_size);
  var_buffer_size = CS_MAX(var_buffer_size, 128);
  var_buffer_size = CS_MIN(var_buffer_size, output_size);

  BFT_MALLOC(var_buffer, var_buffer_size, double);

  /* Export field */
  /*--------------*/

  if (location == FVM_WRITER_PER_NODE) {

    _export_field_n(mesh,
                    helper,
                    writer,
                    writer->bases[base_index - 1],
                    field_label,
                    sol_index,
                    dimension,
                    cgns_datatype,
                    interlace,
                    n_parent_lists,
                    parent_num_shift,
                    datatype,
                    field_values,
                    var_buffer_size,
                    var_buffer);

  } /* End of export field with node location */

  else if (location == FVM_WRITER_PER_ELEMENT) {

    _export_field_e(export_list,
                    helper,
                    writer,
                    writer->bases[base_index - 1],
                    field_label,
                    sol_index,
                    dimension,
                    cgns_datatype,
                    interlace,
                    n_parent_lists,
                    parent_num_shift,
                    datatype,
                    field_values,
                    var_buffer_size,
                    var_buffer);

  } /* End of export_field on elements */

  else
    bft_error(__FILE__, __LINE__, 0,
              "fvm_to_cgns_export_field(): field location not managed.\n"
              "Associated writer: \"%s\"\n"
              "Associated base: \"%s\"\n"
              "Associated field: \"%s\"\n"
              "Associated location: %i\n",
              writer->name, base_name, field_name, location);

  /* Free buffers and helper structures */
  /*------------------------------------*/

  BFT_FREE(var_buffer);

  helper = fvm_writer_field_helper_destroy(helper);

  BFT_FREE(export_list);

  if (rank == 0)
    BFT_FREE(field_label);

}

/*----------------------------------------------------------------------------*/

#endif /* defined(HAVE_CGNS) */

#ifdef __cplusplus
}
#endif /* __cplusplus */
