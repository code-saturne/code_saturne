#ifndef __FVM_TO_CATALYST_H__
#define __FVM_TO_CATALYST_H__

#if defined(HAVE_CATALYST)

/*============================================================================
 * Write a nodal representation associated with a mesh and associated
 * variables to Catalyst objects
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
 * Initialize FVM to Catalyst object writer.
 *
 * Options are:
 *   private_comm        use private MPI communicator (default: false)
 *   names=<fmt>         use same naming rules as <fmt> format
 *                       (default: ensight)
 *
 * parameters:
 *   name           <-- base output case name.
 *   options        <-- whitespace separated, lowercase options list
 *   time_dependecy <-- indicates if and how meshes will change with time
 *   comm           <-- associated MPI communicator.
 *
 * returns:
 *   pointer to opaque Catalyst writer structure.
 *----------------------------------------------------------------------------*/

#if defined(HAVE_MPI)

void *
fvm_to_catalyst_init_writer(const char             *name,
                            const char             *path,
                            const char             *options,
                            fvm_writer_time_dep_t   time_dependency,
                            MPI_Comm                comm);

#else

void *
fvm_to_catalyst_init_writer(const char             *name,
                            const char             *path,
                            const char             *options,
                            fvm_writer_time_dep_t   time_dependency);

#endif

/*----------------------------------------------------------------------------
 * Finalize FVM to Catalyst object writer.
 *
 * parameters:
 *   this_writer_p <-- pointer to opaque Catalyst writer structure.
 *
 * returns:
 *   NULL pointer.
 *----------------------------------------------------------------------------*/

void *
fvm_to_catalyst_finalize_writer(void  *this_writer_p);

/*----------------------------------------------------------------------------
 * Associate new time step with a Catalyst geometry.
 *
 * parameters:
 *   this_writer_p <-- pointer to associated writer
 *   time_step     <-- time step number
 *   time_value    <-- time_value number
 *----------------------------------------------------------------------------*/

void
fvm_to_catalyst_set_mesh_time(void    *this_writer_p,
                              int      time_step,
                              double   time_value);

/*----------------------------------------------------------------------------
 * Write nodal mesh to a a Catalyst object
 *
 * parameters:
 *   this_writer_p <-- pointer to associated writer.
 *   mesh          <-- pointer to nodal mesh structure that should be written.
 *----------------------------------------------------------------------------*/

void
fvm_to_catalyst_export_nodal(void               *this_writer_p,
                             const fvm_nodal_t  *mesh);

/*----------------------------------------------------------------------------
 * Write data to the MultiBlock aimed to be CoProcessed
 *
 * parameters:
 *   this_writer_p    <-- pointer to associated writer
 *   ugrid            <-- pointer to the vtkUnstructuredGrid
 *----------------------------------------------------------------------------*/

void
fvm_to_catalyst_export_field(void                  *this_writer_p,
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
fvm_to_catalyst_flush(void  *this_writer_p);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* defined(HAVE_CATALYST) */

#endif /* __FVM_TO_CATALYST_H__ */
