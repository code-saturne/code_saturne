#ifndef __FVM_TO_CGNS_H__
#define __FVM_TO_CGNS_H__

#if defined(HAVE_CGNS)

/*============================================================================
 * Write a nodal representation associated with a mesh and associated
 * variables to CGNS files
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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
 * Returns number of library version strings associated with the CGNS format.
 *
 * returns:
 *   number of library version strings associated with the CGNS format.
 *----------------------------------------------------------------------------*/

int
fvm_to_cgns_n_version_strings(void);

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
                           int compile_time_version);

/*----------------------------------------------------------------------------
 * Initialize FVM to CGNS file writer.
 *
 * Options are:
 *   discard_polygons    do not output polygons or related values
 *   discard_polyhedra   do not output polyhedra or related values
 *   divide_polygons     tesselate polygons with triangles
 *   adf                 use ADF file type
 *   hdf5                use HDF5 file type (default if available)
 *   links               split output to separate files using links
 *
 * As CGNS does not handle polyhedral elements in a simple manner,
 * polyhedra are automatically tesselated with tetrahedra and pyramids
 * (adding a vertex near each polyhedron's center) unless discarded.
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
                        MPI_Comm                comm);

#else

void *
fvm_to_cgns_init_writer(const char             *name,
                        const char             *path,
                        const char             *options,
                        fvm_writer_time_dep_t   time_dependency);

#endif

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
fvm_to_cgns_finalize_writer(void  *this_writer_p);

/*----------------------------------------------------------------------------
 * Associate new time step with a CGNS geometry.
 *
 * parameters:
 *   this_writer_p <-- pointer to associated writer
 *   time_step     <-- time step number
 *   time_value    <-- time_value number
 *----------------------------------------------------------------------------*/

void
fvm_to_cgns_set_mesh_time(void     *this_writer_p,
                          int       time_step,
                          double    time_value);

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
fvm_to_cgns_needs_tesselation(void               *this_writer_p,
                              const fvm_nodal_t  *mesh,
                              fvm_element_t       element_type);

/*----------------------------------------------------------------------------
 * Write nodal mesh to a CGNS file
 *
 * parameters:
 *   this_writer_p <-- pointer to associated writer.
 *   mesh          <-- pointer to nodal mesh structure that should be written.
 *----------------------------------------------------------------------------*/

void
fvm_to_cgns_export_nodal(void               *this_writer_p,
                         const fvm_nodal_t  *mesh);

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
                         const void       *const field_values[]);

/*----------------------------------------------------------------------------
 * Flush files associated with a given writer.
 *
 * parameters:
 *   this_writer_p    <-- pointer to associated writer
 *----------------------------------------------------------------------------*/

void
fvm_to_cgns_flush(void  *this_writer_p);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* HAVE_CGNS */

#endif /* __FVM_TO_CGNS_H__ */
