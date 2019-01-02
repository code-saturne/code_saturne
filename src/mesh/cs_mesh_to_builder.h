#ifndef __CS_MESH_TO_BUILDER_H__
#define __CS_MESH_TO_BUILDER_H__

/*============================================================================
 * Define cs_mesh_builder_t fields from cs_mesh_t fields.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

#include "fvm_group.h"
#include "fvm_selector.h"
#include "fvm_periodicity.h"

#include "cs_base.h"

#include "cs_io.h"
#include "cs_mesh.h"
#include "cs_mesh_builder.h"
#include "cs_part_to_block.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Transfer mesh to mesh builder structure.
 *
 * As the dataflow is very similar, but may be done array-by array to minimize
 * memory overhead, this function also handles a part of the output
 * to file needed to save a mesh file.
 *
 * parameters:
 *   mesh     <-> pointer to mesh structure
 *   mb       <-> pointer to mesh builder structure
 *   transfer <-- if true, data is transferred from mesh to builder;
 *                if false, builder fields are only used as a temporary
 *                arrays.
 *   pp_out   <-> optional output file, or NULL
 *----------------------------------------------------------------------------*/

void
cs_mesh_to_builder(cs_mesh_t          *mesh,
                   cs_mesh_builder_t  *mb,
                   bool                transfer,
                   cs_io_t            *pp_out);

/*----------------------------------------------------------------------------
 * Transfer mesh partitioning info to mesh builder structure.
 *
 * parameters:
 *   mesh <-- pointer to mesh structure
 *   mb   <-> pointer to mesh builder structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_to_builder_partition(const cs_mesh_t    *mesh,
                             cs_mesh_builder_t  *mb);

/*----------------------------------------------------------------------------
 * Reconstruct periodic faces info from mesh to builder.
 *
 * parameters:
 *   mesh <-- pointer to mesh structure
 *   mb   <-> pointer to mesh builder structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_to_builder_perio_faces(const cs_mesh_t    *mesh,
                               cs_mesh_builder_t  *mb);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MESH_TO_BUILDER_H__ */
