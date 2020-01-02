#ifndef __CS_JOIN_PERIO_H__
#define __CS_JOIN_PERIO_H__

/*============================================================================
 * Structure and function headers handling with periodicity for joining
 * operations
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
#include "cs_mesh.h"
#include "cs_join_set.h"
#include "cs_join_util.h"
#include "cs_join_mesh.h"

/*---------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local Macro definitions
 *===========================================================================*/

/*=============================================================================
 * Global variables
 *===========================================================================*/

/*============================================================================
 *  Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Check if periodic joining operations are queued.
 *
 * Fortran Interface:
 *
 * SUBROUTINE TSTJPE
 * *****************
 *
 * INTEGER        iperio    : <-> : do we have periodicity ?
 * INTEGER        iperot    : <-> : do we have periodicity of rotation ?
 *----------------------------------------------------------------------------*/

void CS_PROCF(tstjpe, tstjpe)
(
 cs_int_t    *iperio,
 cs_int_t    *iperot
);

/*=============================================================================
 * Public function prototypes
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Define a translational periodicity
 *
 * parameters:
 *   sel_criteria  <-- boundary face selection criteria
 *   fraction      <-- value of the fraction parameter
 *   plane         <-- value of the plane parameter
 *   verbosity     <-- level of verbosity required
 *   visualization <-- level of visualization required
 *   trans         <-- translation vector
 *
 * returns:
 *   number (1 to n) associated with new periodicity
 *----------------------------------------------------------------------------*/

int
cs_join_perio_add_translation(const char    *sel_criteria,
                              double         fraction,
                              double         plane,
                              int            verbosity,
                              int            visualization,
                              const double   trans[3]);

/*----------------------------------------------------------------------------
 * Define a rotational periodicity
 *
 * parameters:
 *   sel_criteria  <-- boundary face selection criteria
 *   fraction      <-- value of the fraction parameter
 *   plane         <-- value of the plane parameter
 *   verbosity     <-- level of verbosity required
 *   visualization <-- level of visualization required
 *   theta         <-- rotation angle (in degrees)
 *   axis          <-- axis vector
 *   invariant     <-- invariant point coordinates
 *
 * returns:
 *   joining number (1 to n) associated with new periodicity
 *----------------------------------------------------------------------------*/

int
cs_join_perio_add_rotation(const char    *sel_criteria,
                           double         fraction,
                           double         plane,
                           int            verbosity,
                           int            visualization,
                           double         theta,
                           const double   axis[3],
                           const double   invariant[3]);

/*----------------------------------------------------------------------------
 * Define a periodicity using a matrix
 *
 * parameters:
 *   sel_criteria  <-- boundary face selection criteria
 *   fraction      <-- value of the fraction parameter
 *   plane         <-- value of the plane parameter
 *   verbosity     <-- level of verbosity required
 *   visualization <-- level of visualization required
 *   matrix        <-- transformation matrix
 *
 * returns:
 *   joining number (1 to n) associated with new periodicity
 *----------------------------------------------------------------------------*/

int
cs_join_perio_add_mixed(const char    *sel_criteria,
                        double         fraction,
                        double         plane,
                        int            verbosity,
                        int            visualization,
                        double         matrix[3][4]);

/*----------------------------------------------------------------------------
 * Add periodicity information to mesh and create or update mesh builder
 * for a new periodic joining.
 *
 * parameters:
 *   this_join <-- high level join structure
 *   mesh      <-> pointer to a cs_mesh_t structure
 *   builder   <-> pointer to a cs_mesh_builder_t structure pointer
 *---------------------------------------------------------------------------*/

void
cs_join_perio_init(cs_join_t           *this_join,
                   cs_mesh_t           *mesh,
                   cs_mesh_builder_t  **builder);

/*----------------------------------------------------------------------------
 * Duplicate and apply transformation to the selected faces and also to
 * their related vertices. Modify compact_face_gnum to take into account
 * new periodic faces and create a periodic vertex couple list.
 *
 * parameters:
 *   this_join <-- high level join structure
 *   jmesh     <-> local join mesh struct. to duplicate and transform
 *   mesh      <-- pointer to a cs_mesh_t structure
 *---------------------------------------------------------------------------*/

void
cs_join_perio_apply(cs_join_t          *this_join,
                    cs_join_mesh_t     *jmesh,
                    const cs_mesh_t    *mesh);

/*----------------------------------------------------------------------------
 * Duplicate and apply transformation to the selected faces and also to
 * their related vertices. Modify compact_face_gnum to take into account
 * new periodic faces and create a periodic vertex couple list.
 *
 * parameters:
 *   this_join          <-- pointer to a high level join structure
 *   jmesh              <-> local join mesh struct. to duplicate and transform
 *   mesh               <-- pointer to a cs_mesh_t structure
 *   p_work_jmesh       <-> distributed join mesh struct. on which operations
 *                          take place
 *   p_work_edges       <-> join edges struct. related to work_jmesh
 *   init_max_vtx_gnum  <-- initial max. global numbering for vertices
 *   n_g_new_vertices   <-- global number of vertices created during the
 *                          intersection of edges
 *---------------------------------------------------------------------------*/

void
cs_join_perio_merge_back(cs_join_t          *this_join,
                         cs_join_mesh_t     *jmesh,
                         const cs_mesh_t    *mesh,
                         cs_join_mesh_t    **p_work_jmesh,
                         cs_join_edges_t   **p_work_edges,
                         cs_gnum_t           init_max_vtx_gnum,
                         cs_gnum_t           n_g_new_vertices);

/*----------------------------------------------------------------------------
 * Duplicate and apply transformation to the selected faces and also to
 * their related vertices. Update jmesh structure.
 * Define a new n2o_hist.
 *
 * parameters:
 *   this_join  <-- pointer to a high level join structure
 *   jmesh      <-> local join mesh struct. to duplicate and transform
 *   mesh       <-- pointer to a cs_mesh_t structure
 *   builder    <-- pointer to a cs_mesh_builder_t structure
 *   o2n_hist   <-- old global face -> new local face numbering
 *   p_n2o_hist <-- new global face -> old local face numbering
 *---------------------------------------------------------------------------*/

void
cs_join_perio_split_back(cs_join_t          *this_join,
                         cs_join_mesh_t     *jmesh,
                         cs_mesh_t          *mesh,
                         cs_mesh_builder_t  *builder,
                         cs_join_gset_t     *o2n_hist,
                         cs_join_gset_t    **p_n2o_hist);

/*----------------------------------------------------------------------------
 * Define a list of coupled faces by periodicty in global numbering.
 *
 * For parallel runs:
 *  - remove isolated periodic faces in the mesh definition
 *  - define a consistent face connectivity in order to prepare the building
 *    of periodic vertex couples
 *
 * parameters:
 *   param        <-- set of parameters for the joining operation
 *   n_ii_faces   <-- initial local number of interior faces
 *   face_type    <-- type of faces in join mesh (interior or border ...)
 *   jmesh        <-- pointer to a cs_join_mesh_t structure
 *   mesh         <-> pointer to a cs_mesh_t structure
 *   mesh_builder <-> pointer to a cs_mesh_t structure
 *---------------------------------------------------------------------------*/

void
cs_join_perio_split_update(cs_join_param_t             param,
                           cs_lnum_t                   n_ii_faces,
                           const cs_join_face_type_t   face_type[],
                           const cs_join_mesh_t       *jmesh,
                           cs_mesh_t                  *mesh,
                           cs_mesh_builder_t          *mesh_builder);

/*----------------------------------------------------------------------------
 * Use periodic face couples in cs_glob_join_perio_builder to define
 * cs_glob_mesh_builder
 * Free all elements which can be freed.
 * Transfer data to cs_glob_mesh and cs_glob_mesh_builder.
 *---------------------------------------------------------------------------*/

void
cs_join_perio_transfer_builder(void);

/*---------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_JOIN_PERIO_H__ */
