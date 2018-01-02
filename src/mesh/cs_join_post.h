#ifndef __CS_JOIN_POST_H__
#define __CS_JOIN_POST_H__

/*============================================================================
 * Subroutines to manage post-treatment for joining operation
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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_nodal.h"
#include "fvm_writer.h"

#include "cs_join_mesh.h"
#include "cs_join_util.h"
#include "cs_mesh.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro and type definitions
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create and define a writer to make post-treatment for the joining operation
 *---------------------------------------------------------------------------*/

void
cs_join_post_init(void);

/*----------------------------------------------------------------------------
 * Post-treatment of a cs_join_mesh_t structure.
 *
 * parameters:
 *  mesh_name <-- name of the mesh for the post-processing
 *  mesh      <-- pointer to a cs_join_mesh_t structure to post
 *---------------------------------------------------------------------------*/

void
cs_join_post_mesh(const char            *mesh_name,
                  const cs_join_mesh_t  *join_mesh);

/*----------------------------------------------------------------------------
 * Post-process a subset of faces of a cs_join_mesh_t structure.
 *
 * parameters:
 *   mesh_name        <-- name of the sub-set mesh
 *   mesh             <-- pointer to the parent cs_join_mesh_t structure
 *   n_selected_faces <-- number of selected faces (size of the sub-set)
 *   selected_faces   <-- list of local number in parent mesh
 *---------------------------------------------------------------------------*/

void
cs_join_post_faces_subset(const char            *mesh_name,
                          const cs_join_mesh_t  *parent_mesh,
                          cs_lnum_t              n_select_faces,
                          const cs_lnum_t        selected_faces[]);

/*----------------------------------------------------------------------------
 * Post-process mesh after the update following the fusion operation.
 *
 * parameters:
 *  join_param    <--  set of parameters for the joining operation
 *  join_select   <--  list of all implied entities in the joining operation
 *---------------------------------------------------------------------------*/

void
cs_join_post_after_merge(cs_join_param_t          join_param,
                         const cs_join_select_t  *join_select);


/*----------------------------------------------------------------------------
 * Post-process mesh after the update following the split operation.
 *
 * parameters:
 *  n_old_i_faces   <--  initial number of interior faces
 *  n_old_b_faces   <--  initial number of border faces
 *  n_g_new_b_faces <--  global number of new border faces
 *  n_select_faces  <--  number of selected faces
 *  mesh            <--  pointer to a cs_mesh_t structure
 *  join_param      <--  set of parameters for the joining operation
 *---------------------------------------------------------------------------*/

void
cs_join_post_after_split(cs_lnum_t         n_old_i_faces,
                         cs_lnum_t         n_old_b_faces,
                         cs_gnum_t         n_g_new_b_faces,
                         cs_lnum_t         n_select_faces,
                         const cs_mesh_t  *mesh,
                         cs_join_param_t   join_param);

/*----------------------------------------------------------------------------
 * Post-process mesh after the update following the split operation.
 *
 * parameters:
 *   n_i_clean_faces <-- number of interior faces cleaned
 *   i_clean_faces   <-> list of interior face numbers (ordered on exit)
 *   n_b_clean_faces <-- number of border faces cleaned
 *   b_clean_faces   <-> list of border face numbers (ordered on exit)
 *   param           <-- set of parameters for the joining operation
 *---------------------------------------------------------------------------*/

void
cs_join_post_cleaned_faces(cs_lnum_t        n_i_clean_faces,
                           cs_lnum_t        i_clean_faces[],
                           cs_lnum_t        n_b_clean_faces,
                           cs_lnum_t        b_clean_faces[],
                           cs_join_param_t  param);

/*----------------------------------------------------------------------------
 * Post and dump a cs_join_mesh_t according to the verbosity level.
 *
 * parameters:
 *   basename <--  generic name for the mesh to post
 *   mesh     <--  fvm_join_mesh_t structure to post
 *   param    <--  fvm_join_param_t structure
 *---------------------------------------------------------------------------*/

void
cs_join_post_dump_mesh(const char            *basename,
                       const cs_join_mesh_t  *mesh,
                       cs_join_param_t        param);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_JOIN_POST_H__ */
