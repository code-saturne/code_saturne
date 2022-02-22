#ifndef __FVM_WRITER_H__
#define __FVM_WRITER_H__

/*============================================================================
 * Handle export of mesh and fields.
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

#include "cs_timer.h"

#include "fvm_defs.h"
#include "fvm_nodal.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Element types
 *----------------------------------------------------------------------------*/

typedef enum {

  FVM_WRITER_FIXED_MESH,         /* Mesh definitions do not change with time */
  FVM_WRITER_TRANSIENT_COORDS,   /* Vertex coordinates may change with time */
  FVM_WRITER_TRANSIENT_CONNECT   /* Mesh connectivity may change with time */

} fvm_writer_time_dep_t;

/*----------------------------------------------------------------------------
 * Variable definition type
 *----------------------------------------------------------------------------*/

typedef enum {

  FVM_WRITER_PER_NODE,           /* Variable values per node */
  FVM_WRITER_PER_ELEMENT,        /* Variable values per element */
  FVM_WRITER_PER_PARTICLE        /* Variable values per particle */

} fvm_writer_var_loc_t;

/*----------------------------------------------------------------------------
 * Opaque structure defining a writer definition
 *----------------------------------------------------------------------------*/

typedef struct _fvm_writer_t fvm_writer_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/* Names of time dependency enumeration values */

extern const char  *fvm_writer_time_dep_name[];

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Find the format matching a name,
 *
 * parameters:
 *   format_name <-- name of desired format
 *
 * returns:
 *   index of the format matching the given name, or -1 if none matches.
 *----------------------------------------------------------------------------*/

int
fvm_writer_get_format_id(const char  *format_name);

/*----------------------------------------------------------------------------
 * Returns number of known formats.
 *----------------------------------------------------------------------------*/

int
fvm_writer_n_formats(void);

/*----------------------------------------------------------------------------
 * Returns name of a known format.
 *
 * parameters:
 *   format_index <-- index of format in known format list (0 to n-1)
 *
 * returns:
 *   pointer to constant string containing the format's name
 *----------------------------------------------------------------------------*/

const char *
fvm_writer_format_name(int format_index);

/*----------------------------------------------------------------------------
 * Returns availability of a known format.
 *
 * parameters:
 *   format_index <-- index of format in known format list (0 to n-1)
 *
 * returns:
 *   1 if the format is available, 0 otherwise.
 *----------------------------------------------------------------------------*/

int
fvm_writer_format_available(int format_index);

/*----------------------------------------------------------------------------
 * Returns number of library version strings associated with a given format.
 *
 * For writers requiring an external library, the first associated
 * version string should correspond to that library, with possible
 * additional version strings for its dependencies.
 *
 * For writers only requiring standard libraries (libc, MPI, MPI-IO),
 * this function should return 0.
 *
 * parameters:
 *   format_index <-- index of format in known format list (0 to n-1)
 *
 * returns:
 *   number of library version strings associated with a given format.
 *----------------------------------------------------------------------------*/

int
fvm_writer_n_version_strings(int format_index);

/*----------------------------------------------------------------------------
 * Returns a library version string associated with a given format.
 *
 * We must have string_index < fvm_writer_n_version_strings(format_index).
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
 *   format_index <-- index of format in known format list (0 to n-1)
 *   string_index <-- index in format's version string list (0 to n-1)
 *   compile_time <-- 0 by default, 1 if we want the compile-time version
 *                    string, if different from the run-time version.
 *
 * returns:
 *   pointer to constant string containing the library's version.
 *----------------------------------------------------------------------------*/

const char *
fvm_writer_version_string(int format_index,
                          int string_index,
                          int compile_time_version);

/*----------------------------------------------------------------------------
 * Remove a given option from a format options list if present.
 *
 * The possible separators are also transformed to whitespace and merged.
 *
 * parameters:
 *   format_options  <-> options for the selected format (case-independent,
 *                       whitespace or comma separated list)
 *   exclude_option  <-- option to be excluded
 *----------------------------------------------------------------------------*/

void
fvm_writer_filter_option(char        *format_options,
                         const char  *exclude_option);

/*----------------------------------------------------------------------------
 * Initialize FVM mesh and field output writer.
 *
 * Allowed options depend on what is applicable to a given format. Those
 * not relevant to a given writer are ignored. Possible options include:
 *   text                output text files (EnSight)
 *   binary              output binary files (EnSight, default)
 *   big_endian          force binary files to big-endian (EnSight)
 *   adf                 use ADF file type (CGNS)
 *   hdf5                use HDF5 file type (CGNS, default if available)
 *   discard_polygons    do not output polygons or related values
 *   discard_polyhedra   do not output polyhedra or related values
 *   divide_polygons     tesselate polygons with triangles
 *   divide_polyhedra    tesselate polyhedra with tetrahedra and pyramids
 *                       (adding a vertex near each polyhedron's center)
 *   separate_meshes     use a different writer for each mesh
 *
 * parameters:
 *   name            <-- base name of output
 *   path            <-- optional directory name for output
 *                       (directory automatically created if necessary)
 *   format_name     <-- name of selected format (case-independent)
 *   format_options  <-- options for the selected format (case-independent,
 *                       whitespace or comma separated list)
 *   time_dependency <-- indicates if and how meshes will change with time
 *
 * returns:
 *   pointer to mesh and field output writer
 *----------------------------------------------------------------------------*/

fvm_writer_t *
fvm_writer_init(const char             *name,
                const char             *path,
                const char             *format_name,
                const char             *format_options,
                fvm_writer_time_dep_t   time_dependency);

/*----------------------------------------------------------------------------
 * Finalize FVM mesh and field output writer.
 *
 * parameters:
 *   this_writer <-- pointer to mesh and field output writer
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

fvm_writer_t *
fvm_writer_finalize(fvm_writer_t  *this_writer);

/*----------------------------------------------------------------------------
 * Return a writer's name.
 *
 * parameters:
 *   this_writer <-- pointer to mesh and field output writer
 *
 * returns:
 *   pointer to base name of output associated with the writer
 *----------------------------------------------------------------------------*/

const char *
fvm_writer_get_name(const fvm_writer_t  *this_writer);

/*----------------------------------------------------------------------------
 * Return a writer's associated format name.
 *
 * parameters:
 *   this_writer <-- pointer to mesh and field output writer
 *
 * returns:
 *   pointer to output format name associated with the writer
 *----------------------------------------------------------------------------*/

const char *
fvm_writer_get_format(const fvm_writer_t  *this_writer);

/*----------------------------------------------------------------------------
 * Return a writer's associated format options.
 *
 * parameters:
 *   this_writer <-- pointer to mesh and field output writer
 *
 * returns:
 *   pointer to output format options associated with the writer
 *----------------------------------------------------------------------------*/

const char *
fvm_writer_get_options(const fvm_writer_t  *this_writer);

/*----------------------------------------------------------------------------
 * Return a writer's associated output directory.
 *
 * parameters:
 *   this_writer <-- pointer to mesh and field output writer
 *
 * returns:
 *   pointer to output format options associated with the writer
 *----------------------------------------------------------------------------*/

const char *
fvm_writer_get_path(const fvm_writer_t  *this_writer);

/*----------------------------------------------------------------------------
 * Return geometry time dependency status of a writer.
 *
 * parameters:
 *   this_writer <-- pointer to mesh and field output writer
 *
 * returns:
 *   time dependency status
 *----------------------------------------------------------------------------*/

fvm_writer_time_dep_t
fvm_writer_get_time_dep(const fvm_writer_t  *this_writer);

/*----------------------------------------------------------------------------
 * Associate new time step with a mesh.
 *
 * parameters:
 *   this_writer_p <-- pointer to associated writer
 *   time_step     <-- time step number
 *   time_value    <-- time_value number
 *----------------------------------------------------------------------------*/

void
fvm_writer_set_mesh_time(fvm_writer_t  *this_writer,
                         int            time_step,
                         double         time_value);

/*----------------------------------------------------------------------------
 * Query if elements of a given type will need to be tesselated
 * for use of a nodal mesh with an output writer.
 *
 * This function should be called before any fvm_writer_export_...()
 *
 * parameters:
 *   this_writer  <-- pointer to mesh and field output writer
 *   mesh         <-- pointer to nodal mesh
 *   element_type <-- type of element
 *
 * returns:
 *   0 if no tesselation is necessary, 1 if tesselation is necessary.
 *----------------------------------------------------------------------------*/

int
fvm_writer_needs_tesselation(fvm_writer_t       *this_writer,
                             const fvm_nodal_t  *mesh,
                             fvm_element_t       element_type);

/*----------------------------------------------------------------------------
 * Export FVM nodal mesh.
 *
 * parameters:
 *   this_writer <-- pointer to mesh and field output writer
 *   mesh        <-- pointer to nodal mesh
 *----------------------------------------------------------------------------*/

void
fvm_writer_export_nodal(fvm_writer_t       *this_writer,
                        const fvm_nodal_t  *mesh);

/*----------------------------------------------------------------------------
 * Export field associated with a nodal mesh.
 *
 * Assigning a negative value to the time step indicates a time-independent
 * field (in which case the time_value argument is unused).
 *
 * parameters:
 *   this_writer      <-- pointer to mesh and field output writer
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
fvm_writer_export_field(fvm_writer_t                 *this_writer,
                        const fvm_nodal_t            *mesh,
                        const char                   *name,
                        fvm_writer_var_loc_t          location,
                        int                           dimension,
                        cs_interlace_t                interlace,
                        int                           n_parent_lists,
                        const cs_lnum_t               parent_num_shift[],
                        cs_datatype_t                 datatype,
                        int                           time_step,
                        double                        time_value,
                        const void             *const field_values[]);

/*----------------------------------------------------------------------------
 * Flush files associated with a given writer.
 *
 * parameters:
 *   this_writer      <-- pointer to mesh and field output writer
 *----------------------------------------------------------------------------*/

void
fvm_writer_flush(fvm_writer_t  *this_writer);

/*----------------------------------------------------------------------------
 * Return accumulated times associated with output for a given writer.
 *
 * parameters:
 *   this_writer <-- pointer to mesh and field output writer
 *   mesh_time   --> Meshes output time (or NULL)
 *   field_time  --> Fields output time (or NULL)
 *   flush_time  --> remaining (applying output) time (or NULL)
 *----------------------------------------------------------------------------*/

void
fvm_writer_get_times(fvm_writer_t        *this_writer,
                     cs_timer_counter_t  *mesh_time,
                     cs_timer_counter_t  *field_time,
                     cs_timer_counter_t  *flush_time);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __FVM_WRITER_H__ */
