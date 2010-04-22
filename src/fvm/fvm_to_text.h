#ifndef __FVM_TO_TEXT_H__
#define __FVM_TO_TEXT_H__

/*============================================================================
 * Write a nodal representation associated with a mesh to file
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2004-2006  EDF

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_defs.h"
#include "fvm_nodal.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

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
 * Initialize FVM to text file writer.
 *
 * parameters:
 *   name    <-- base output case name.
 *   options <-- whitespace separated, lowercase options list
 *   comm    <-- associated MPI communicator.
 *
 * returns:
 *   pointer to opaque text file writer structure.
 *----------------------------------------------------------------------------*/

#if defined(HAVE_MPI)

void *
fvm_to_text_init_writer(const char                   *const name,
                        const char                   *const path,
                        const char                   *const options,
                        const fvm_writer_time_dep_t         time_dependency,
                        const MPI_Comm                      comm);

#else

void *
fvm_to_text_init_writer(const char                   *const name,
                        const char                   *const path,
                        const char                   *const options,
                        const fvm_writer_time_dep_t         time_dependency);

#endif

/*----------------------------------------------------------------------------
 * Finalize FVM to text file writer.
 *
 * parameters:
 *   this_writer_p <-- pointer to opaque text file writer structure.
 *
 * returns:
 *   NULL pointer.
 *----------------------------------------------------------------------------*/

void *
fvm_to_text_finalize_writer(void  *this_writer_p);

/*----------------------------------------------------------------------------
 * Write nodal mesh to a text file
 *
 * parameters:
 *   this_writer_p <-- pointer to associated writer.
 *   mesh          <-- pointer to nodal mesh structure that should be written.
 *----------------------------------------------------------------------------*/

void
fvm_to_text_export_nodal(void               *const this_writer_p,
                         const fvm_nodal_t  *const mesh);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __FVM_TO_TEXT_H__ */
