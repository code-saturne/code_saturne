#ifndef __ECS_COMM_H__
#define __ECS_COMM_H__

/*============================================================================
 * Base functions for writing Kernel I/O files.
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

#include "ecs_def.h"

BEGIN_C_DECLS

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Opaque structure to handle file output */

typedef struct _ecs_comm_t ecs_comm_t;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize a Kernel I/O file writer.
 *
 * returns:
 *   initialized Kernel I/O file writer
 *----------------------------------------------------------------------------*/

ecs_comm_t *
ecs_comm_initialize(const char  *file_name);

/*----------------------------------------------------------------------------
 * Close writer.
 *
 * arguments:
 *   comm <-- pointer to writer structure pointer
 *----------------------------------------------------------------------------*/

void
ecs_comm_finalize(ecs_comm_t **comm);

/*----------------------------------------------------------------------------
 * Write a section to the Kernel I/O file.
 *
 * Grid locations and possibly indexes may be assigned to a section by
 * specifying a location id; the first time a given location id appears in
 * the file is considered a declaration. In the same manner, an index id
 * may be specified. Values of zero indicate no location or index base
 * is used. It is up to the calling code to ensure that total number of
 * values, location size, and number of values per location are consistent,
 * as this my be important for code reading the file.
 *
 * arguments:
 *   name              <-- section name
 *   location_id       <-- id of associated location
 *   index_id          <-- id of associated index
 *   n_location_values <-- number of values per location
 *   embed             <-- embed values in header
 *   values            <-- values to write
 *   value_type        <-- type of value to write
 *   comm              <-- Kernel I/O file output structure
 *----------------------------------------------------------------------------*/

void
ecs_comm_write_section(const char  *name,
                       size_t       n_values,
                       size_t       location_id,
                       size_t       index_id,
                       size_t       n_location_values,
                       bool         embed,
                       const void  *values,
                       ecs_type_t   value_type,
                       ecs_comm_t  *comm);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __ECS_COMM_H__ */
