#ifndef __FVM_TO_CATALYST2_H__
#define __FVM_TO_CATALYST2_H__

#if defined(HAVE_CATALYST2)

/*============================================================================
 * Write a nodal representation associated with a mesh and associated
 * variables to Catalyst2 objects
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

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

#include "base/cs_defs.h"

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
 * Returns number of library version strings associated with Catalyst-2.
 *
 * The first associated version string should corresponds to Catalyst,
 * The second to the Conduit library.
 *
 * returns:
 *   number of library version strings associated with Catalyst output.
 *----------------------------------------------------------------------------*/

int
fvm_to_catalyst2_n_version_strings(void);

/*----------------------------------------------------------------------------
 * Returns a library version string associated with the Catalyst-2 output.
 *
 * The first associated version string should correspond to Catalyst info,
 * The second to conduit info.
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
fvm_to_catalyst2_version_string(int string_index,
                                int compile_time_version);

/*----------------------------------------------------------------------------
 * Initialize FVM to Catalyst-2 object writer.
 *
 * Options are:
 *   private_comm        use private MPI communicator (default: false)
 *   names=<fmt>         use same naming rules as <fmt> format
 *                       (default: ensight)
 *   input_name=<name>   define input name (default: writer name)
 *
 * parameters:
 *   name           <-- base output case name.
 *   options        <-- whitespace separated, lowercase options list
 *   time_dependecy <-- indicates if and how meshes will change with time
 *   comm           <-- associated MPI communicator.
 *
 * returns:
 *   pointer to opaque writer structure.
 *----------------------------------------------------------------------------*/

#if defined(HAVE_MPI)
void *
fvm_to_catalyst2_init_writer(const char             *name,
                             const char             *path,
                             const char             *options,
                             fvm_writer_time_dep_t   time_dependency,
                             MPI_Comm                comm);
#else
void *
fvm_to_catalyst2_init_writer(const char             *name,
                             const char             *path,
                             const char             *options,
                             fvm_writer_time_dep_t   time_dependency);

#endif

/*----------------------------------------------------------------------------
 * Finalize FVM to Catalyst-2 object writer.
 *
 * parameters:
 *   this_writer_p <-- pointer to opaque writer structure.
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

void *
fvm_to_catalyst2_finalize_writer(void  *this_writer_p);

/*----------------------------------------------------------------------------
 * Associate new time step with a Catalyst-2 geometry.
 *
 * parameters:
 *   this_writer_p <-- pointer to associated writer
 *   time_step     <-- time step number
 *   time_value    <-- time_value number
 *----------------------------------------------------------------------------*/

void
fvm_to_catalyst2_set_mesh_time(void    *this_writer_p,
                               int      time_step,
                               double   time_value);

/*----------------------------------------------------------------------------
 * Write nodal mesh to a Catalyst-2 object
 *
 * parameters:
 *   this_writer_p <-- pointer to associated writer.
 *   mesh          <-- pointer to nodal mesh structure that should be written.
 *----------------------------------------------------------------------------*/

void
fvm_to_catalyst2_export_nodal(void               *this_writer_p,
                              const fvm_nodal_t  *mesh);

/*----------------------------------------------------------------------------
 * Write field associated with a nodal mesh to to a Catalyst-2 object.
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
fvm_to_catalyst2_export_field(void                  *this_writer_p,
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
                              const void      *const field_values[]);

/*----------------------------------------------------------------------------
 * Flush files associated with a given writer.
 *
 * In this case, the effective call to coprocessing is done.
 *
 * parameters:
 *   this_writer_p    <-- pointer to associated writer
 *----------------------------------------------------------------------------*/

void
fvm_to_catalyst2_flush(void  *this_writer_p);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* defined(HAVE_CATALYST2) */

#endif /* __FVM_TO_CATALYST2_H__ */
