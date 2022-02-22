#ifndef __CS_MESH_FROM_BUILDER_H__
#define __CS_MESH_FROM_BUILDER_H__

/*============================================================================
 * Define cs_mesh_t fields from cs_mesh_builder_t fields.
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

#include "fvm_group.h"
#include "fvm_selector.h"
#include "fvm_periodicity.h"

#include "cs_base.h"

#include "cs_mesh.h"
#include "cs_mesh_builder.h"

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
 * Transfer mesh builder to mesh structure.
 *
 * parameters:
 *   mesh         <-- pointer to mesh structure
 *   mesh_builder <-- pointer to mesh builder structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_from_builder(cs_mesh_t          *mesh,
                     cs_mesh_builder_t  *mesh_builder);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MESH_FROM_BUILDER_H__ */
