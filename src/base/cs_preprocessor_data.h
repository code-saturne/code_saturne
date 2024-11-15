#ifndef __CS_PRE_PROCESSOR_DATA_H__
#define __CS_PRE_PROCESSOR_DATA_H__

/*============================================================================
 * Exchange of data between code_saturne Kernel and Preprocessor
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_io_num.h"

#include "cs_base.h"
#include "cs_mesh.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! Mesh preprocessing operations in case of restart. */

typedef enum {

  CS_PREPROCESSOR_DATA_RESTART_NONE,        /*!< do not use restart mesh
                                              for computation, or restart
                                              mesh not present */
  CS_PREPROCESSOR_DATA_RESTART_AND_MODIFY,  /*!< read restart mesh and enable
                                               further preprocessing */
  CS_PREPROCESSOR_DATA_RESTART_ONLY,        /*!< use restart mesh as-is,
                                              with no additional
                                              preprocessing (default if
                                              restart is present) */

} cs_preprocessor_data_restart_mode_t;

/*============================================================================
 *  Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Pass information relative to mesh metadata to the Fortran API
 *
 * Fortran Interface:
 *
 * subroutine ledevi(iperio, iperot)
 * *****************
 *
 * integer          iperio      : <-- : Periodicity indicator
 *----------------------------------------------------------------------------*/

void
CS_PROCF(ledevi, LEDEVI)(int  *iperio);

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
 *   transf_matrix   <-- coordinate transformation matrix (or null)
 *----------------------------------------------------------------------------*/

void
cs_preprocessor_data_add_file(const char     *file_name,
                              size_t          n_group_renames,
                              const char    **group_rename,
                              const double    transf_matrix[3][4]);

/*----------------------------------------------------------------------------
 * Check for periodicity information in mesh meta-data.
 *
 * returns:
 *   0 if no periodicity is present in mesh input,
 *   1 for translation periodicity only,
 *   2 for rotation or mixed periodicity
 *----------------------------------------------------------------------------*/

int
cs_preprocessor_check_perio(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return restart behavior for preprocessing.
 *
 * \return  preprocessing mode in case of restart.
 */
/*----------------------------------------------------------------------------*/

cs_preprocessor_data_restart_mode_t
cs_preprocessor_data_get_restart_mode(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define restart behavior in case of restart.
 *
 * If no restart/mesh_input.csm (or restart/mesh_input) file is found,
 * CS_PREPROCESSOR_DATA_RESTART_NONE will be used.
 *
 * \param[in]  mode  chosen preprocessing mode on restart
 */
/*----------------------------------------------------------------------------*/

void
cs_preprocessor_data_set_restart_mode(cs_preprocessor_data_restart_mode_t  mode);

/*----------------------------------------------------------------------------
 * Read mesh meta-data.
 *
 * parameters:
 *   mesh             <-- pointer to mesh structure
 *   mesh_builder     <-- pointer to mesh builder structure
 *   ignore_cartesian <-- option to ignore cartesian blocks
 *----------------------------------------------------------------------------*/

void
cs_preprocessor_data_read_headers(cs_mesh_t          *mesh,
                                  cs_mesh_builder_t  *mesh_builder,
                                  bool                ignore_cartesian);

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
 *   mesh             <-- pointer to mesh structure
 *   mesh_builder     <-- pointer to mesh builder structure
 *   ignore_cartesian <-- option to ignore cartesian blocks
 *----------------------------------------------------------------------------*/

void
cs_preprocessor_data_read_mesh(cs_mesh_t          *mesh,
                               cs_mesh_builder_t  *mesh_builder,
                               bool                ignore_cartesian);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PRE_PROCESSOR_DATA_H__ */

