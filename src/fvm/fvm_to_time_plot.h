#ifndef __FVM_TO_TIME_PLOT_H__
#define __FVM_TO_TIME_PLOT_H__

/*============================================================================
 * Write a nodal representation associated with a mesh and associated
 * variables to time plot (whitepace-separated or csv) files
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
 * Initialize FVM to time plot file writer.
 *
 * Options are:
 *   csv                 output CSV (comma-separated-values) files
 *   dat                 output dat (space-separated) files
 *   use_iteration       use time step id instead of time value for
 *                       first column
 *   flush_wtime=<wt>    flush output file every 'wt' seconds
 *   n_buf_steps=<n>     write output to file every 'n' output steps
 *
 * parameters:
 *   name           <-- base output case name.
 *   options        <-- whitespace separated, lowercase options list
 *   time_dependecy <-- indicates if and how meshes will change with time
 *   comm           <-- associated MPI communicator.
 *
 * returns:
 *   pointer to opaque time plot writer structure.
 *----------------------------------------------------------------------------*/

#if defined(HAVE_MPI)

void *
fvm_to_time_plot_init_writer(const char             *name,
                             const char             *path,
                             const char             *options,
                             fvm_writer_time_dep_t   time_dependency,
                             MPI_Comm                comm);

#else

void *
fvm_to_time_plot_init_writer(const char             *name,
                             const char             *path,
                             const char             *options,
                             fvm_writer_time_dep_t  time_dependency);

#endif

/*----------------------------------------------------------------------------
 * Finalize FVM to time plot file writer.
 *
 * parameters:
 *   writer <-- pointer to opaque time plot writer structure.
 *
 * returns:
 *   NULL pointer.
 *----------------------------------------------------------------------------*/

void *
fvm_to_time_plot_finalize_writer(void  *writer);

/*----------------------------------------------------------------------------
 * Associate new time step with a time plot geometry.
 *
 * parameters:
 *   writer     <-- pointer to associated writer
 *   time_step  <-- time step number
 *   time_value <-- time_value number
 *----------------------------------------------------------------------------*/

void
fvm_to_time_plot_set_mesh_time(void          *writer,
                               const int      time_step,
                               const double   time_value);

/*----------------------------------------------------------------------------
 * Write nodal mesh to a time plot file
 *
 * parameters:
 *   writer <-- pointer to associated writer.
 *   mesh   <-- pointer to nodal mesh structure that should be written.
 *----------------------------------------------------------------------------*/

void
fvm_to_time_plot_export_nodal(void               *writer,
                              const fvm_nodal_t  *mesh);

/*----------------------------------------------------------------------------
 * Write field associated with a nodal mesh to a time plot file.
 *
 * Assigning a negative value to the time step indicates a time-independent
 * field (in which case the time_value argument is unused).
 *
 * parameters:
 *   writer           <-- pointer to associated writer
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
fvm_to_time_plot_export_field(void                  *writer,
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

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __FVM_TO_TIME_PLOT_H__ */
