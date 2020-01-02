#ifndef __CS_JOIN_UPDATE_H__
#define __CS_JOIN_UPDATE_H__

/*============================================================================
 * Structure and function headers handling with mesh update during
 * the joining operation
 *===========================================================================*/

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

/*----------------------------------------------------------------------------
 *  Local headers
 *---------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_join_mesh.h"
#include "cs_join_set.h"
#include "cs_join_util.h"
#include "cs_mesh.h"

/*---------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *===========================================================================*/

/*============================================================================
 * Type definition
 *===========================================================================*/

/*=============================================================================
 * Global variables
 *===========================================================================*/

/*============================================================================
 *  Public function header for Fortran API
 *===========================================================================*/

/*=============================================================================
 * Public function prototypes
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Update mesh structure (vertices + faces) after the fusion step.
 *
 * parameters:
 *   join_param   <--  set of parameters for the joining operation
 *   join_select  <--  list of all implied entities in the joining op.
 *   o2n_vtx_gnum <->  in : array on blocks on the new global vertex
 *                     out: local array on the new global vertex
 *   join_mesh    <->  pointer to the local cs_join_mesh_t structure
 *   mesh         <->  pointer of pointer to cs_mesh_t structure
 *---------------------------------------------------------------------------*/

void
cs_join_update_mesh_after_merge(cs_join_param_t    join_param,
                                cs_join_select_t  *join_select,
                                cs_gnum_t          o2n_vtx_gnum[],
                                cs_join_mesh_t    *join_mesh,
                                cs_mesh_t         *mesh);

/*----------------------------------------------------------------------------
 * Update mesh structure (vertices + faces) after the face split step.
 *
 * parameters:
 *   join_param      <-- set of parameters for the joining operation
 *   join_select     <-- list of all implied entities in the joining op.
 *   o2n_face_hist   <-- relation between faces before/after the joining
 *   join_mesh       <-> pointer to the local cs_join_mesh_t structure
 *   mesh            <-> pointer of pointer to cs_mesh_t structure
 *   mesh_builder    <-> pointer of pointer to cs_mesh__builder_t structure
 *---------------------------------------------------------------------------*/

void
cs_join_update_mesh_after_split(cs_join_param_t          join_param,
                                const cs_join_select_t  *join_select,
                                const cs_join_gset_t    *o2n_face_hist,
                                cs_join_mesh_t          *join_mesh,
                                cs_mesh_t               *mesh,
                                cs_mesh_builder_t       *mesh_builder);

/*----------------------------------------------------------------------------
 * Clean a cs_mesh_t struct.
 * Delete redundant and empty edge definitions.
 *
 * parameters:
 *   para    <-- set of parameters for the joining operation
 *   mesh    <-> pointer to a cs_mesh_t structure
 *---------------------------------------------------------------------------*/

void
cs_join_update_mesh_clean(cs_join_param_t   param,
                          cs_mesh_t        *mesh);

/*---------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_JOIN_UPDATE_H__ */
