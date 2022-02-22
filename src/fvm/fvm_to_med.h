#ifndef __FVM_TO_MED_H__
#define __FVM_TO_MED_H__

#if defined(HAVE_MED)

/*============================================================================
 * Write a nodal representation associated with a mesh and associated
 * variables to MED files
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_defs.h"
#include "fvm_nodal.h"
#include "fvm_writer.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
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
fvm_to_med_n_version_strings(void);

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
                          int compile_time_version);

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
fvm_to_med_init_writer(const char                   *const name,
                       const char                   *const path,
                       const char                   *const options,
                       const fvm_writer_time_dep_t         time_dependency,
                       const MPI_Comm                      comm);

#else

void *
fvm_to_med_init_writer(const char                   *const name,
                       const char                   *const path,
                       const char                   *const options,
                       const fvm_writer_time_dep_t         time_dependency);

#endif

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
fvm_to_med_finalize_writer(void  *this_writer_p);

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
                             fvm_element_t       element_type);

/*----------------------------------------------------------------------------
 * Associate new time step with a MED mesh.
 *
 * parameters:
 *   this_writer <-- pointer to associated writer
 *   time_step   <-- time step number
 *   time_value  <-- time_value number
 *----------------------------------------------------------------------------*/

void
fvm_to_med_set_mesh_time(void          *const this_writer,
                         const int            time_step,
                         const double         time_value);

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
                     const fvm_nodal_t  *mesh);

/*----------------------------------------------------------------------------
 * Write nodal mesh to a MED file
 *
 * parameters:
 *   this_writer_p <-- pointer to associated writer.
 *   mesh          <-- pointer to nodal mesh structure that should be written.
 *----------------------------------------------------------------------------*/

void
fvm_to_med_export_nodal(void               *const this_writer_p,
                        const fvm_nodal_t  *const mesh);

/*----------------------------------------------------------------------------
 * Write field associated with a nodal mesh to a MED file.
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
fvm_to_med_export_field(void                      *const this_writer,
                        const fvm_nodal_t         *const mesh,
                        const char                *const name,
                        const fvm_writer_var_loc_t       location,
                        const int                        dimension,
                        const cs_interlace_t             interlace,
                        const int                        n_parent_lists,
                        const cs_lnum_t                  parent_num_shift[],
                        const cs_datatype_t              datatype,
                        const int                        time_step,
                        const double                     time_value,
                        const void                *const field_values[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __FVM_TO_MED_H__ */

#endif /* HAVE_MED */
