/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 2008-2010 EDF S.A., France
 *
 *     contact: saturne-support@edf.fr
 *
 *     The Code_Saturne Kernel is free software; you can redistribute it
 *     and/or modify it under the terms of the GNU General Public License
 *     as published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version.
 *
 *     The Code_Saturne Kernel is distributed in the hope that it will be
 *     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with the Code_Saturne Kernel; if not, write to the
 *     Free Software Foundation, Inc.,
 *     51 Franklin St, Fifth Floor,
 *     Boston, MA  02110-1301  USA
 *
 *============================================================================*/

#ifndef __CS_JOIN_PERIO_H__
#define __CS_JOIN_PERIO_H__

/*============================================================================
 * Structure and function headers handling with periodicity for joining
 * operations
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * FVM library headers
 *---------------------------------------------------------------------------*/

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

extern int  cs_glob_n_join_perio;   /* Number of periodicity defined through
                                       a joining operation */

/*============================================================================
 *  Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Get the number of periodic transformations already defined
 *
 * Fortran Interface:
 *
 * SUBROUTINE NUMPER
 * *****************
 *
 * INTEGER   numper      : --> : number of periodicities  op. already defined
 *----------------------------------------------------------------------------*/

void CS_PROCF(numper, NUMPER)
(
 cs_int_t    *numper
);

/*----------------------------------------------------------------------------
 * Define a translation
 *
 * Fortran Interface:
 *
 * SUBROUTINE DEFPT1
 * *****************
 *
 * INTEGER        numper    : <-- : number related to the periodic op.
 * CHARACTER*     criteria  : <-- : boundary face selection criteria
 * REAL           fraction  : <-- : parameter for merging vertices
 * REAL           plane     : <-- : parameter for splitting faces
 * INTEGER        verbosity : <-- : verbosity level
 * REAL           tx        : <-- : X coordinate of the translation vector
 * REAL           ty        : <-- : Y coordinate of the translation vector
 * REAL           tz        : <-- : Z coordinate of the translation vector
 * INTEGER        crit_len  : <-- : length of criteria
 *----------------------------------------------------------------------------*/

void CS_PROCF(defpt1, DEFPT1)
(
 cs_int_t    *numper,
 const char  *criteria,
 cs_real_t   *fraction,
 cs_real_t   *plane,
 cs_int_t    *verbosity,
 cs_real_t   *tx,
 cs_real_t   *ty,
 cs_real_t   *tz,
 cs_int_t    *crit_len
 CS_ARGF_SUPP_CHAINE
);

/*----------------------------------------------------------------------------
 * Define a rotation
 *
 * Fortran Interface:
 *
 * SUBROUTINE DEFPR1
 * *****************
 *
 * INTEGER        numper    : <-- : number related to the periodic op.
 * CHARACTER*     criteria : <-- : boundary face selection criteria
 * REAL           fraction  : <-- : parameter for merging vertices
 * REAL           plane     : <-- : parameter for splitting faces
 * INTEGER        verbosity : <-- : verbosity level
 * REAL           ax       : <-- : X coordinate of the rotation axis
 * REAL           ay       : <-- : Y coordinate of the rotation axis
 * REAL           az       : <-- : Z coordinate of the rotation axis
 * REAL           theta    : <-- : angle of the rotation (radian)
 * REAL           ix       : <-- : X coordinate of the invariant point
 * REAL           iy       : <-- : Y coordinate of the invariant point
 * REAL           iz       : <-- : Z coordinate of the invariant point
 * INTEGER        crit_len : <-- : length of criteria string
 *----------------------------------------------------------------------------*/

void CS_PROCF(defpr1, DEFPR1)
(
 cs_int_t    *numper,
 const char  *criteria,
 cs_real_t   *fraction,
 cs_real_t   *plane,
 cs_int_t    *verbosity,
 cs_real_t   *ax,
 cs_real_t   *ay,
 cs_real_t   *az,
 cs_real_t   *theta,
 cs_real_t   *ix,
 cs_real_t   *iy,
 cs_real_t   *iz,
 cs_int_t    *crit_len
 CS_ARGF_SUPP_CHAINE
);

/*----------------------------------------------------------------------------
 * Define a general transformation through a homogeneous matrix (4x4)
 *     _               _
 *    | r11 r12 r13 tx  |  t(x,y,z) : translation vector
 *    | r21 r22 r23 ty  |  r(i,j)   : rotation matrix
 *    | r31 r32 r33 tz  |
 *    |_  0   0   0  1 _|
 *
 * Fortran Interface:
 *
 * SUBROUTINE DEFPG1
 * *****************
 *
 * INTEGER        numper    : <-- : number related to the periodic op.
 * CHARACTER*     criteria  : <-- : boundary face selection criteria
 * REAL           fraction  : <-- : parameter for merging vertices
 * REAL           plane     : <-- : parameter for splitting faces
 * INTEGER        verbosity : <-- : verbosity level
 * REAL           r11       : <-- : coef. (1,1) of the homogeneous matrix
 * REAL           r12       : <-- : coef. (1,2) of the homogeneous matrix
 * REAL           r13       : <-- : coef. (1,3) of the homogeneous matrix
 * REAL           tx        : <-- : coef. (1,4) of the homogeneous matrix
 * REAL           r21       : <-- : coef. (2,1) of the homogeneous matrix
 * REAL           r22       : <-- : coef. (2,2) of the homogeneous matrix
 * REAL           r23       : <-- : coef. (2,3) of the homogeneous matrix
 * REAL           ty        : <-- : coef. (2,4) of the homogeneous matrix
 * REAL           r31       : <-- : coef. (3,1) of the homogeneous matrix
 * REAL           r32       : <-- : coef. (3,2) of the homogeneous matrix
 * REAL           r33       : <-- : coef. (3,3) of the homogeneous matrix
 * REAL           tz        : <-- : coef. (3,4) of the homogeneous matrix
 * INTEGER        crit_len  : <-- : length of criteria string
 *----------------------------------------------------------------------------*/

void CS_PROCF(defpg1, DEFPG1)
(
 cs_int_t    *numper,
 const char  *criteria,
 cs_real_t   *fraction,
 cs_real_t   *plane,
 cs_int_t    *verbosity,
 cs_real_t   *r11,
 cs_real_t   *r12,
 cs_real_t   *r13,
 cs_real_t   *tx,
 cs_real_t   *r21,
 cs_real_t   *r22,
 cs_real_t   *r23,
 cs_real_t   *ty,
 cs_real_t   *r31,
 cs_real_t   *r32,
 cs_real_t   *r33,
 cs_real_t   *tz,
 cs_int_t    *crit_len
 CS_ARGF_SUPP_CHAINE
);

/*----------------------------------------------------------------------------
 * Set advanced parameters for the joining algorithm in case of periodicity
 *
 * Fortran Interface:
 *
 * SUBROUTINE SETAPP
 * *****************
 *
 * INTEGER      perio_num         : <-- : perio number
 * REAL         mtf               : <-- : merge tolerance coefficient
 * REAL         pmf               : <-- : pre-merge factor
 * INTEGER      tcm               : <-- : tolerance computation mode
 * INTEGER      icm               : <-- : intersection computation mode
 * INTEGER      maxbrk            : <-- : max number of equiv. breaks
 * INTEGER      max_sub_faces     : <-- : max. possible number of sub-faces
 *                                        by splitting a selected face
 * INTEGER      tml               : <-- : tree max level
 * INTEGER      tmb               : <-- : tree max boxes
 * REAL         tmr               : <-- : tree max ratio
 *----------------------------------------------------------------------------*/

void CS_PROCF(setapp, SETAPP)
(
 cs_int_t    *perio_num,
 cs_real_t   *mtf,
 cs_real_t   *pmf,
 cs_int_t    *tcm,
 cs_int_t    *icm,
 cs_int_t    *maxbrk,
 cs_int_t    *max_sub_faces,
 cs_int_t    *tml,
 cs_int_t    *tmb,
 cs_real_t   *tmr
);

/*=============================================================================
 * Public function prototypes
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Define a translational periodicity
 *
 * parameters:
 *   perio_num    <-- number related to the periodicity
 *   sel_criteria <-- boundary face selection criteria
 *   fraction     <-- value of the fraction parameter
 *   plane        <-- value of the plane parameter
 *   verbosity    <-- level of verbosity required
 *   trans        <-- translation vector
 *----------------------------------------------------------------------------*/

void
cs_join_perio_add_translation(int            perio_num,
                              const char    *sel_criteria,
                              double         fraction,
                              double         plane,
                              int            verbosity,
                              const double   trans[3]);

/*----------------------------------------------------------------------------
 * Define a rotational periodicity
 *
 * parameters:
 *   perio_num    <-- number related to the periodicity
 *   sel_criteria <-- boundary face selection criteria
 *   fraction     <-- value of the fraction parameter
 *   plane        <-- value of the plane parameter
 *   verbosity    <-- level of verbosity required
 *   theta        <-- rotation angle (in degrees)
 *   axis         <-- axis vector
 *   invariant    <-- invariant point coordinates
 *----------------------------------------------------------------------------*/

void
cs_join_perio_add_rotation(int            perio_num,
                           const char    *sel_criteria,
                           double         fraction,
                           double         plane,
                           int            verbosity,
                           double         theta,
                           const double   axis[3],
                           const double   invariant[3]);

/*----------------------------------------------------------------------------
 * Define a periodicity using a matrix
 *
 * parameters:
 *   perio_num    <-- number related to the periodicity
 *   sel_criteria <-- boundary face selection criteria
 *   fraction     <-- value of the fraction parameter
 *   plane        <-- value of the plane parameter
 *   verbosity    <-- level of verbosity required
 *   matrix       <-- transformation matrix
 *----------------------------------------------------------------------------*/

void
cs_join_perio_add_mixed(int            perio_num,
                        const char    *sel_criteria,
                        double         fraction,
                        double         plane,
                        int            verbosity,
                        double         matrix[3][4]);

/*----------------------------------------------------------------------------
 * Duplicate and apply transformation to the selected faces and also to
 * their related vertices. Modify compact_face_gnum to take into account
 * new periodic faces and create a periodic vertex couple list.
 *
 * parameters:
 *   this_join    <--  high level join structure
 *   jmesh        <->  local join mesh struct. to duplicate and transform
 *   mesh         <--  pointer to a cs_mesh_t struct.
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
                         cs_join_mesh_t    **p_work_jmesh,
                         cs_join_edges_t   **p_work_edges,
                         fvm_gnum_t          init_max_vtx_gnum,
                         fvm_gnum_t          n_g_new_vertices);

/*----------------------------------------------------------------------------
 * Duplicate and apply transformation to the selected faces and also to
 * their related vertices. Update jmesh structure.
 * Define a new n2o_hist.
 *
 * parameters:
 *   this_join     <-- pointer to a high level join structure
 *   jmesh         <-> local join mesh struct. to duplicate and transform
 *   mesh          <-- pointer to a cs_mesh_t structure
 *   o2n_hist      <-- old global face -> new local face numbering
 *   p_n2o_hist    <-- new global face -> old local face numbering
 *---------------------------------------------------------------------------*/

void
cs_join_perio_split_back(cs_join_t          *this_join,
                         cs_join_mesh_t     *jmesh,
                         cs_mesh_t          *mesh,
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
 *
 * parameters:
 *   param        <-- set of parameters for the joining operation
 *   n_ii_faces   <-- initial local number of interior faces
 *   face_type    <-- type of faces in join mesh (interior or border ...)
 *   jmesh        <-- pointer on a cs_join_mesh_t struct.
 *   mesh         <-> pointer on a cs_mesh_t struct.
 *---------------------------------------------------------------------------*/

void
cs_join_perio_split_update(cs_join_param_t             param,
                           cs_int_t                    n_ii_faces,
                           const cs_join_face_type_t   face_type[],
                           const cs_join_mesh_t       *jmesh,
                           cs_mesh_t                  *mesh);

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
