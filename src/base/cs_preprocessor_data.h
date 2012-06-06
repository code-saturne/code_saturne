#ifndef __CS_PRE_PROCESSOR_DATA_H__
#define __CS_PRE_PROCESSOR_DATA_H__

/*============================================================================
 * Exchange of data between Code_Saturne Kernel and Preprocessor
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2012 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include "fvm_io_num.h"

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_mesh.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 *  Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Receive messages from the pre-processor about the dimensions of mesh
 * parameters
 *
 * Fortran Interface:
 *
 * subroutine ledevi(ndim   , nfml  , nprfml, iperio, iperot)
 * *****************
 *
 * integer          ndim        : --> : Spacial dimension (3)
 * integer          nfml        : <-- : Number of families
 * integer          nprfml      : <-- : Number of properties per family
 * integer          iperio      : <-- : Periodicity indicator
 * integer          iperot      : <-- : Number of rotation periodicities
 *----------------------------------------------------------------------------*/

void
CS_PROCF(ledevi, LEDEVI)(const cs_int_t   *ndim,
                         cs_int_t         *nfml,
                         cs_int_t         *nprfml,
                         cs_int_t         *iperio,
                         cs_int_t         *iperot);

/*============================================================================
 *  Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define input mesh file to read.
 *
 * If this function is never called, the default file is read.
 * The first time this function is called,  this default is overriden by the
 * defined file, and all subsequent calls define additional meshes to read.
 *
 * parameters:
 *   file_name       <-- name of file to read
 *   n_group_renames <-- number of groups to rename
 *   group_rename    <-- old (group_rename[i*2]) to new (group_rename[i*2 + 1])
 *                       group names array (size: n_group_renames*2)
 *   transf_matrix   <-- coordinate transformation matrix (or NULL)
 *----------------------------------------------------------------------------*/

void
cs_preprocessor_data_add_file(const char     *file_name,
                              size_t          n_group_renames,
                              const char    **group_rename,
                              const double    transf_matrix[3][4]);

/*----------------------------------------------------------------------------
 * Read pre-processor mesh data and finalize input.
 *
 * At this stage, ghost cells are not generated yet, so the interior
 * face connectivity is not complete near parallel domain or periodic
 * boundaries. Also, isolated faces, if present, are considered to be
 * boundary faces, as they may participate in future mesh joining
 * operations. Their matching cell number will be set to -1.
 * Remaining isolated faces should be removed before completing
 * the mesh structure.
 *
 * parameters:
 *   mesh         <-- pointer to mesh structure
 *   mesh_builder <-- pointer to mesh builder structure
 *----------------------------------------------------------------------------*/

void
cs_preprocessor_data_read_mesh(cs_mesh_t          *mesh,
                               cs_mesh_builder_t  *mesh_builder);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PRE_PROCESSOR_DATA_H__ */

