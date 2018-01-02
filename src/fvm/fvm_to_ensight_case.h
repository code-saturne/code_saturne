#ifndef __FVM_TO_ENSIGHT_CASE_H__
#define __FVM_TO_ENSIGHT_CASE_H__

/*============================================================================
 * Manage case files associated with the EnSight Gold writer
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

/* Opaque structure to manage case file */

typedef struct _fvm_to_ensight_case_t fvm_to_ensight_case_t;

/* Geometry or variable file info */

typedef struct {

  const char  * name;    /* Pointer to file name */
  bool          queried; /* Indicates if this file name has already been
                            returned by "fvm_to_ensight_case_get_..._file()"
                            (so we can decide to create or append to it) */

} fvm_to_ensight_case_file_info_t;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create a new case file structure.
 *
 * parameters:
 *   name            <-- case name
 *   dir_prefix      <-- associated local or absolute directory name
 *   time_dependency <-- indicates if and how meshes will change with time
 *
 * returns:
 *   pointer to new case file structure
 *----------------------------------------------------------------------------*/

fvm_to_ensight_case_t *
fvm_to_ensight_case_create(const char                   *const name,
                           const char                   *const dir_prefix,
                           const fvm_writer_time_dep_t         time_dependency);

/*----------------------------------------------------------------------------
 * Destroy a case file structure.
 *
 * parameters:
 *   this_case  <-- case structure
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

fvm_to_ensight_case_t *
fvm_to_ensight_case_destroy(fvm_to_ensight_case_t  *this_case);

/*----------------------------------------------------------------------------
 * Return time dependency status of an EnSight geometry.
 *
 * parameters:
 *   this_case  <-- case structure
 *
 * returns:
 *   time dependency status
 *----------------------------------------------------------------------------*/

fvm_writer_time_dep_t
fvm_to_ensight_case_get_time_dep(fvm_to_ensight_case_t  *this_case);

/*----------------------------------------------------------------------------
 * Associate new time step with an EnSight geometry.
 *
 * parameters:
 *   this_case  <-- case structure
 *   time_step  <-- time step number
 *   time_value <-- time_value number
 *
 * returns:
 *   0 if no time was added, 1 if a new time was added
 *----------------------------------------------------------------------------*/

int
fvm_to_ensight_case_set_geom_time(fvm_to_ensight_case_t  *const this_case,
                                  const int                     time_step,
                                  const double                  time_value);

/*----------------------------------------------------------------------------
 * Return current file name and "queried" indicator associated with an
 * EnSight geometry.
 *
 * The "queried" flag in the info structure is set to "false" the first
 * time this function returns a given file name, and to "true" all other
 * times.
 *
 * parameters:
 *   this_case  <-- case structure
 *
 * returns:
 *   Info structure for geometry file
 *----------------------------------------------------------------------------*/

fvm_to_ensight_case_file_info_t
fvm_to_ensight_case_get_geom_file(fvm_to_ensight_case_t  *const this_case);

/*----------------------------------------------------------------------------
 * Associate a part name with a case and return its number.
 * If the part was already associated, zero is returned.
 *
 * parameters:
 *   this_case  <-- case structure
 *   part_name  <-- part name
 *
 * returns:
 *   part number in case, or 0 if part already associated
 *----------------------------------------------------------------------------*/

int
fvm_to_ensight_case_add_part(fvm_to_ensight_case_t  *const this_case,
                             const char             *const part_name);

/*----------------------------------------------------------------------------
 * Return the part number associated with a given part name, or 0
 *
 * parameters:
 *   this_case  <-- case structure
 *   part_name  <-- part name
 *
 * returns:
 *   part number in case, or 0 if part name is not associated with this case
 *----------------------------------------------------------------------------*/

int
fvm_to_ensight_case_get_part_num(fvm_to_ensight_case_t  *const this_case,
                                 const char             *const part_name);

/*----------------------------------------------------------------------------
 * Return current file name and "queried" indicator associated with an
 * EnSight variable.
 *
 * The "queried" flag in the info structure is set to "false" the first
 * time this function returns a given file name, and to "true" all other
 * times.
 *
 * if the corresponding variable or physical time are not present in the
 * structure, the necessary elements are added.
 *
 * parameters:
 *   this_case  <-> pointer to structure that should be updated
 *   name       <-- variable name
 *   dimension  <-- variable dimension (0: constant, 1: scalar, 3: vector,
 *                  6: symmetrical tensor, 9: asymmetrical tensor)
 *   location   <-- variable definition location (nodes, elements, or particles)
 *   time_step  <-- number of time step to add
 *   time_value <-- associated time value
 *
 * returns:
 *   Info structure for file associated with the variable
 *----------------------------------------------------------------------------*/

fvm_to_ensight_case_file_info_t
fvm_to_ensight_case_get_var_file(fvm_to_ensight_case_t       *const this_case,
                                 const char                  *const name,
                                 const int                          dimension,
                                 const fvm_writer_var_loc_t         location,
                                 const int                          time_step,
                                 const double                       time_value);

/*----------------------------------------------------------------------------
 * Write an EnSight Gold case file.
 *
 * parameters:
 *   this_case  <-- case structure
 *   rank       <-- calling rank in case of parallelism
 *----------------------------------------------------------------------------*/

void
fvm_to_ensight_case_write_case(fvm_to_ensight_case_t  *this_case,
                               int                     rank);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __FVM_TO_ENSIGHT_CASE_H__ */
