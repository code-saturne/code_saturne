#ifndef __CS_POST_H__
#define __CS_POST_H__

/*============================================================================
 * Post-processing management
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2011 EDF S.A.

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
 * BFT library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <fvm_nodal.h>
#include <fvm_writer.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*
 * Output type masks
 */

#define CS_POST_VOLUME               (1 << 0)  /* postprocess volume */
#define CS_POST_BOUNDARY_NR          (1 << 1)  /* postprocess bundary
                                                  without reconstruction */

/*============================================================================
 * Local type definitions
 *============================================================================*/

/* Datatype enumeration */

typedef enum {
  CS_POST_TYPE_cs_int_t,
  CS_POST_TYPE_cs_real_t,
  CS_POST_TYPE_int,
  CS_POST_TYPE_float,
  CS_POST_TYPE_double
} cs_post_type_t;

/* Function pointer associated with a specific post-processing variables
   output: such functions are registered using the cs_post_add_time_dep_var(),
   and all registered functions are automatically called by PSTVAR. */

typedef void
(cs_post_time_dep_var_t) (int        instance_id,
                          int        nt_cur_abs,
                          cs_real_t  t_cur_abs);

/*=============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public Fortran function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Configure the post-processing output so that a mesh displacement field
 * may be output automatically for meshes based on the global volume mesh/
 *
 * Fortran interface:
 *
 * subroutine pstdfm
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstdfm, PSTDFM)
(
 void
);

/*----------------------------------------------------------------------------
 * Update the "active" or "inactive" flag for writers based on the current
 * time step and their default output frequency.
 *
 * Fortran interface:
 *
 * subroutine pstntc (ntmabs, ntcabs, ttcabs)
 * *****************
 *
 * integer          ntmabs      : <-- : maximum time step number
 * integer          ntcabs      : <-- : current time step number
 * double precision ttcabs      : <-- : absolute time at the current time step
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstntc, PSTNTC)
(
 const cs_int_t  *ntmabs,
 const cs_int_t  *ntcabs,
 const cs_real_t *ttcabs
);

/*----------------------------------------------------------------------------
 * Force the "active" or "inactive" flag for a specific writer or for all
 * writers for the current time step.
 *
 * Fortran interface:
 *
 * subroutine pstact (numwri, indact)
 * *****************
 *
 * integer          numwri      : <-- : writer number, or 0 for all writers
 * integer          indact      : <-- : 0 to deactivate, 1 to activate
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstact, PSTACT)
(
 const cs_int_t  *numwri,
 const cs_int_t  *indact
);

/*----------------------------------------------------------------------------
 * Output post-processing meshes using associated writers.
 *
 * Fortran interface:
 *
 * subroutine pstema (ntcabs, ttcabs)
 * *****************
 *
 * integer          ntcabs      : <-- : current time step number
 * double precision ttcabs      : <-- : current physical time
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstema, PSTEMA)
(
 const cs_int_t   *ntcabs,
 const cs_real_t  *ttcabs
);

/*----------------------------------------------------------------------------
 * Loop on post-processing meshes to output variables
 *
 * Fortran interface:
 *
 * subroutine pstvar
 * *****************
 *                  ( ntcabs,
 *                    nvar,   nscal,  nvlsta, nvisbr,
 *                    ttcabs,
 *                    dt,     rtpa,   rtp,    propce, propfa, propfb,
 *                    coefa,  coefb,
 *                    statce, stativ, statfb,
 *                    ra)
 *
 * integer          ntcabs      : --> : current time step number
 * integer          nvar        : <-- : number of variables
 * integer          nscal       : <-- : number of scalars
 * integer          nvlsta      : <-- : number of statistical variables (lagr)
 * integer          nvisbr      : <-- : number of boundary stat. variables (lagr)
 * double precision ttcabs      : <-- : current physical time
 * double precision dt          : <-- : local time step
 * double precision rtpa        : <-- : cell variables at previous time step
 * double precision rtp         : <-- : cell variables
 * double precision propce      : <-- : cell physical properties
 * double precision propfa      : <-- : interior face physical properties
 * double precision propfb      : <-- : boundary face physical properties
 * double precision coefa       : <-- : boundary conditions array
 * double precision coefb       : <-- : boundary conditions array
 * double precision statce      : <-- : cell statistics (lagrangian)
 * double precision stativ      : <-- : cell variance statistics (lagrangian)
 * double precision statfb      : <-- : boundary face statistics (lagrangian)
 * double precision ra          : <-- : ra floating-point array
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstvar, PSTVAR)
(
 const cs_int_t   *ntcabs,
 const cs_int_t   *nvar,
 const cs_int_t   *nscal,
 const cs_int_t   *nvlsta,
 const cs_int_t   *nvisbr,
 const cs_real_t  *ttcabs,
 const cs_real_t   dt[],
 const cs_real_t   rtpa[],
 const cs_real_t   rtp[],
 const cs_real_t   propce[],
 const cs_real_t   propfa[],
 const cs_real_t   propfb[],
 const cs_real_t   coefa[],
 const cs_real_t   coefb[],
 const cs_real_t   statce[],
 const cs_real_t   stativ[],
 const cs_real_t   statfb[],
       cs_real_t   ra[]
);

/*----------------------------------------------------------------------------
 * Post-processing output of a variable defined on cells or faces of a mesh
 * using associated writers.
 *
 * fortran interface; use psteva (see cs_post_f2c.f90)
 *
 * subroutine pstev1 (nummai, nomvar, lnmvar, idimt,  ientla, ivarpr,
 * *****************
 *                    ntcabs, ttcabs, varcel, varfac, varfbr)
 *
 * integer          nummai      : <-- : number of associated output mesh
 * character        nomvar      : <-- : name of associated variable
 * integer          lnmvar      : <-- : variable name length
 * integer          idimt       : <-- : 1 for scalar, 3 for vector
 * integer          ientla      : <-- : if a vector, 1 for interlaced values
 *                              :     : (x1, y1, z1, x2, y2, ..., yn, zn),
 *                              :     : 0 otherwise (x1, x2, ...xn, y1, y2, ...)
 * integer          ivarpr      : <-- : 1 if variable is defined on "parent"
 *                              :     : mesh, 2 if defined on output mesh
 * integer          ntcabs      : <-- : current time step number
 * double precision ttcabs      : <-- : current physical time
 * double precision varcel(*)   : <-- : cell values
 * double precision varfac(*)   : <-- : interior face values
 * double precision varfbo(*)   : <-- : boundary face values
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstev1, PSTEV1)
(
 const cs_int_t   *nummai,
 const char       *nomvar,
 const cs_int_t   *lnmvar,
 const cs_int_t   *idimt,
 const cs_int_t   *ientla,
 const cs_int_t   *ivarpr,
 const cs_int_t   *ntcabs,
 const cs_real_t  *ttcabs,
 const cs_real_t   varcel[],
 const cs_real_t   varfac[],
 const cs_real_t   varfbr[]
 CS_ARGF_SUPP_CHAINE              /*     (possible 'length' arguments added
                                         by many Fortran compilers) */
);

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Configure the post-processing output so that a mesh displacement field
 * may be output automatically for meshes based on the global volume mesh/
 *----------------------------------------------------------------------------*/

void
cs_post_set_deformable(void);

/*----------------------------------------------------------------------------
 * Define a writer; this objects manages a case's name, directory, and format,
 * as well as associated mesh's time dependency, and the default output
 * frequency for associated variables.
 *
 * This function must be called before the time loop. If a writer with a
 * given id is defined multiple times, the last definition supercedes the
 * previous ones.
 *
 * parameters:
 *   writer_id     <-- number of writer to create (< 0 reserved, > 0 for user)
 *   case_name     <-- associated case name
 *   dir_name      <-- associated directory name
 *   fmt_name      <-- associated format name
 *   fmt_opts      <-- associated format options string
 *   time_dep      <-- FVM_WRITER_FIXED_MESH if mesh definitions are fixed,
 *                     FVM_WRITER_TRANSIENT_COORDS if coordinates change,
 *                     FVM_WRITER_TRANSIENT_CONNECT if connectivity changes
 *   output_at_end <-- force output at calculation end if not 0
 *   frequency_n   <-- default output frequency in time-steps, or < 0
 *   frequency_t   <-- default output frequency in seconds, or < 0
 *                     (has priority over frequency_n)
 *----------------------------------------------------------------------------*/

void
cs_post_define_writer(int                     writer_id,
                      const char             *case_name,
                      const char             *dir_name,
                      const char             *fmt_name,
                      const char             *fmt_opts,
                      fvm_writer_time_dep_t   time_dep,
                      bool                    output_at_end,
                      cs_int_t                frequency_n,
                      cs_real_t               frequency_t);

/*----------------------------------------------------------------------------
 * Define a volume post-processing mesh.
 *
 * parameters:
 *   mesh_id        <-- id of mesh to define (< 0 reserved, > 0 for user)
 *   mesh_name      <-- associated mesh name
 *   cell_criteria  <-- selection criteria for cells
 *   add_groups     <-- if true, add group information if present
 *   auto_variables <-- if true, automatic output of main variables
 *   n_writers      <-- number of associated writers
 *   writer_ids     <-- ids of associated writers
 *----------------------------------------------------------------------------*/

void
cs_post_define_volume_mesh(int          mesh_id,
                           const char  *mesh_name,
                           const char  *cell_criteria,
                           bool         add_groups,
                           bool         auto_variables,
                           int          n_writers,
                           const int    writer_ids[]);

/*----------------------------------------------------------------------------
 * Define a volume post-processing mesh using a cell list.

 * The list of cells to extract is sorted upon exit, whether it was sorted
 * upon calling or not.
 *
 * parameters:
 *   mesh_id        <-- id of mesh to define (< 0 reserved, > 0 for user)
 *   mesh_name      <-- associated mesh name
 *   n_cells        <-- number of selected cells
 *   cell_list      <-> list of selected cells (1 to n numbering)
 *   add_groups     <-- if true, add group information if present
 *   auto_variables <-- if true, automatic output of main variables
 *   n_writers      <-- number of associated writers
 *   writer_ids     <-- ids of associated writers
 *----------------------------------------------------------------------------*/

void
cs_post_define_volume_mesh_by_list(int          mesh_id,
                                   const char  *mesh_name,
                                   cs_lnum_t    n_cells,
                                   cs_lnum_t    cell_list[],
                                   bool         add_groups,
                                   bool         auto_variables,
                                   int          n_writers,
                                   const int    writer_ids[]);

/*----------------------------------------------------------------------------
 * Define a surface post-processing mesh.
 *
 * parameters:
 *   mesh_id         <-- id of mesh to define (< 0 reserved, > 0 for user)
 *   mesh_name       <-- associated mesh name
 *   i_face_criteria <-- selection criteria for interior faces
 *   b_face_criteria <-- selection criteria for boundary faces
 *   add_groups      <-- if true, add group information if present
 *   auto_variables  <-- if true, automatic output of main variables
 *   n_writers       <-- number of associated writers
 *   writer_ids      <-- ids of associated writers
 *----------------------------------------------------------------------------*/

void
cs_post_define_surface_mesh(int          mesh_id,
                            const char  *mesh_name,
                            const char  *i_face_criteria,
                            const char  *b_face_criteria,
                            bool         add_groups,
                            bool         auto_variables,
                            int          n_writers,
                            const int    writer_ids[]);

/*----------------------------------------------------------------------------
 * Define a surface post-processing mesh using a face list.
 *
 * Lists of cells or faces to extract are sorted upon exit, whether they
 * were sorted upon calling or not.
 *
 * parameters:
 *   mesh_id        <-- id of mesh to define (< 0 reserved, > 0 for user)
 *   mesh_name      <-- associated mesh name
 *   n_i_faces      <-- number of associated interior faces
 *   n_b_faces      <-- number of associated boundary faces
 *   i_face_list    <-> list of associated interior faces (1 to n numbering)
 *   b_face_list    <-> list of associated boundary faces (1 to n numbering)
 *   add_groups     <-- if true, add group information if present
 *   auto_variables <-- if true, automatic output of main variables
 *   n_writers      <-- number of associated writers
 *   writer_ids     <-- ids of associated writers
 *----------------------------------------------------------------------------*/

void
cs_post_define_surface_mesh_by_list(int          mesh_id,
                                    const char  *mesh_name,
                                    cs_lnum_t    n_i_faces,
                                    cs_lnum_t    n_b_faces,
                                    cs_lnum_t    i_face_list[],
                                    cs_lnum_t    b_face_list[],
                                    bool         add_groups,
                                    bool         auto_variables,
                                    int          n_writers,
                                    const int    writer_ids[]);

/*----------------------------------------------------------------------------
 * Create an alias to a post-processing mesh.
 *
 * An alias allows association of an extra identifier (id) to an
 * existing post-processing mesh, and thus to associate different writers
 * than those associated with the existing mesh. For example, this allows
 * outputting a set of main variables every n1 time steps with one writer,
 * and outputting a specific set of variables every n2 time time steps to
 * another post-processing set using another writer, without the overhead
 * that would be incurred by duplication of the post-processing mesh.
 *
 * An alias is thus treated in all points like its associated mesh;
 * if the definition of either one is modified, that of the other is
 * modified also.
 *
 * It is forbidden to associate an alias to another alias (as there is no
 * identified use for this, and it would make consistency checking more
 * difficult), but multiple aliases may be associated with a given mesh.
 *
 * parameters:
 *   mesh_id         <-- id of mesh to define (< 0 reserved, > 0 for user)
 *   aliased_mesh_id <-- id of aliased mesh
 *   auto_variables  <-- if true, automatic output of main variables
 *   n_writers       <-- number of associated writers
 *   writer_ids      <-- ids of associated writers
 *----------------------------------------------------------------------------*/

void
cs_post_define_alias_mesh(int        mesh_id,
                          int        aliased_mesh_id,
                          bool       auto_variables,
                          int        n_writers,
                          const int  writer_ids[]);

/*----------------------------------------------------------------------------
 * Create a post-processing mesh associated with an existing exportable mesh
 * representation.
 *
 * If the exportable mesh is not intended to be used elsewhere, one can choose
 * to transfer its property to the post-processing mesh, which will then
 * manage its lifecycle based on its own requirements.
 *
 * If the exportable mesh must still be shared, one must be careful to
 * maintain consistency between this mesh and the post-processing output.
 *
 * The mesh in exportable dimension may be of a lower dimension than
 * its parent mesh, if it has been projected. In this case, a
 * dim_shift value of 1 indicates that parent cells are mapped to
 * exportable faces, and faces to edges, while a dim_shift value of 2
 * would indicate that parent cells are mapped to edges.
 * This is important when variables values are exported.
 *
 * parameters:
 *   mesh_id        <-- number of mesh to create (< 0 reserved, > 0 for user)
 *   exp_mesh       <-- mesh in exportable representation (i.e. fvm_nodal_t)
 *   dim_shift      <-- nonzero if exp_mesh has been projected
 *   transfer       <-- if true, ownership of exp_mesh is transferred to
 *                      the post-processing mesh
 *   auto_variables <-- if true, automatic output of main variables
 *   n_writers      <-- number of associated writers
 *   writer_ids     <-- ids of associated writers
 *----------------------------------------------------------------------------*/

void
cs_post_define_existing_mesh(int           mesh_id,
                             fvm_nodal_t  *exp_mesh,
                             int           dim_shift,
                             bool          transfer,
                             bool          auto_variables,
                             int           n_writers,
                             const int     writer_ids[]);

/*----------------------------------------------------------------------------
 * Create a mesh based upon the extraction of edges from an existing mesh.
 *
 * The newly created edges have no link to their parent elements, so
 * no variable referencing parent elements may be output to this mesh,
 * whose main use is to visualize "true" face edges when polygonal faces
 * are subdivided by the writer. In this way, even highly non-convex
 * faces may be visualized correctly if their edges are overlaid on
 * the surface mesh with subdivided polygons.
 *
 * parameters:
 *   mesh_id <-- id of edges mesh to create (< 0 reserved, > 0 for user)
 *   base_mesh_id   <-- id of existing mesh (< 0 reserved, > 0 for user)
 *   n_writers  <-- number of associated writers
 *   writer_ids <-- ids of associated writers
 *----------------------------------------------------------------------------*/

void
cs_post_define_edges_mesh(int        mesh_id,
                          int        base_mesh_id,
                          int        n_writers,
                          const int  writer_ids[]);

/*----------------------------------------------------------------------------
 * Remove a post-processing mesh.
 *
 * No further post-processing output will be allowed on this mesh,
 * so the associated structures may be freed.
 *
 * A post-processing mesh that has been associated with a time-varying
 * writer or that is referenced by an alias may not be removed.
 *
 * parameters:
 *   mesh_id <-- id of mesh to remove
 *----------------------------------------------------------------------------*/

void
cs_post_free_mesh(int  mesh_id);

/*----------------------------------------------------------------------------
 * Check for the existence of a writer of the given id.
 *
 * parameters:
 *   writer_id <-- writer id to check
 *
 * returns:
 *   true if writer with this id exists, false otherwise
 *----------------------------------------------------------------------------*/

bool
cs_post_writer_exists(int  writer_id);

/*----------------------------------------------------------------------------
 * Return a pointer to the FVM library writer associated to a writer_id.
 *
 * parameters:
 *   writer_id <-- associated writer id
 *
 * Returns:
 *  a pointer to a fvm_writer_t structure
 *----------------------------------------------------------------------------*/

fvm_writer_t *
cs_post_get_writer(cs_int_t  writer_id);

/*----------------------------------------------------------------------------
 * Check for the existence of a post-processing mesh of the given id.
 *
 * parameters:
 *   mesh_id <-- mesh id to check
 *
 * returns:
 *   true if mesh with this id exists, false otherwise
 *----------------------------------------------------------------------------*/

bool
cs_post_mesh_exists(int  mesh_id);

/*----------------------------------------------------------------------------
 * Modify an existing post-processing mesh.
 *
 * The lists of cells or faces are redefined, for example to update an
 * extracted mesh based in "interesting" zones.
 *
 * It is not necessary to use this function if a mesh is simply deformed.
 *
 * parameters:
 *   mesh_id     <-- id of mesh to modify (< 0 reserved, > 0 for user)
 *   n_cells     <-- number of associated cells
 *   n_i_faces   <-- number of associated interior faces
 *   n_b_faces   <-- number of associated boundary faces
 *   cell_list   <-> list of associated cells
 *   i_face_list <-> list of associated interior faces
 *   b_face_list <-> list of associated boundary faces
 *
 *----------------------------------------------------------------------------*/

void
cs_post_modify_mesh(int       mesh_id,
                    cs_int_t  n_cells,
                    cs_int_t  n_i_faces,
                    cs_int_t  n_b_faces,
                    cs_int_t  cell_list[],
                    cs_int_t  i_face_list[],
                    cs_int_t  b_face_list[]);

/*----------------------------------------------------------------------------
 * Return the default writer format name
 *
 * Returns:
 *   name of the default writer format
 *----------------------------------------------------------------------------*/

const char *
cs_post_get_default_format(void);

/*----------------------------------------------------------------------------
 * Return the default writer format options
 *
 * Returns:
 *   default writer format options string
 *----------------------------------------------------------------------------*/

const char *
cs_post_get_default_format_options(void);

/*----------------------------------------------------------------------------
 * Return the next "reservable" (i.e. non-user) writer id available.
 *
 * Returns:
 *   the smallest negative integer present, -1
 *----------------------------------------------------------------------------*/

int
cs_post_get_free_writer_id(void);

/*----------------------------------------------------------------------------
 * Return the next "reservable" (i.e. non-user) mesh id available.
 *
 * Returns:
 *   the smallest negative integer present, -1
 *----------------------------------------------------------------------------*/

int
cs_post_get_free_mesh_id(void);

/*----------------------------------------------------------------------------
 * Update "active" or "inactive" flag of writers whose output frequency
 * is a divisor of the current time step number.
 *
 * parameters:
 *   nt_max_abs <-- maximum time step number
 *   nt_cur_abs <-- current time step number
 *   t_cur_abs  <-- absolute time at the current time step
 *----------------------------------------------------------------------------*/

void
cs_post_activate_if_default(int     nt_max_abs,
                            int     nt_cur_abs,
                            double  t_cur_abs);

/*----------------------------------------------------------------------------
 * Force the "active" or "inactive" flag for a specific writer or for all
 * writers for the current time step.
 *
 * parameters:
 *   writer_id <-- writer id, or 0 for all writers
 *   activate  <-- false to deactivate, true to activate
 *----------------------------------------------------------------------------*/

void
cs_post_activate_writer(int   writer_id,
                        bool  activate);

/*----------------------------------------------------------------------------
 * Output post-processing meshes using associated writers.
 *
 * parameters:
 *   nt_cur_abs <-- current time step number
 *   t_cur_abs  <-- current physical time
 *----------------------------------------------------------------------------*/

void
cs_post_write_meshes(int     nt_cur_abs,
                     double  t_cur_abs);

/*----------------------------------------------------------------------------
 * Output a variable defined at cells or faces of a post-processing mesh
 * using associated writers.
 *
 * parameters:
 *   mesh_id     <-- id of associated mesh
 *   var_name    <-- name of variable to output
 *   var_dim     <-- 1 for scalar, 3 for vector
 *   interlace   <-- if a vector, true for interlaced values, false otherwise
 *   use_parent  <-- true if values are defined on "parent" mesh,
 *                   false if values are defined on post-processing mesh
 *   var_type    <-- variable's data type
 *   nt_cur_abs  <-- current time step number
 *   t_cur_abs   <-- current physical time
 *   cel_vals    <-- cell values
 *   i_face_vals <-- interior face values
 *   b_face_vals <-- boundary face values
 *----------------------------------------------------------------------------*/

void
cs_post_write_var(int              mesh_id,
                  const char      *var_name,
                  cs_int_t         var_dim,
                  bool             interlace,
                  bool             use_parent,
                  cs_post_type_t   var_type,
                  cs_int_t         nt_cur_abs,
                  cs_real_t        t_cur_abs,
                  const void      *cel_vals,
                  const void      *i_face_vals,
                  const void      *b_face_vals);

/*----------------------------------------------------------------------------
 * Output a variable defined at vertices of a post-processing mesh using
 * associated writers.
 *
 * parameters:
 *   mesh_id    <-- id of associated mesh
 *   var_name   <-- name of variable to output
 *   var_dim    <-- 1 for scalar, 3 for vector
 *   interlace  <-- if a vector, true for interlaced values, false otherwise
 *   use_parent <-- true if values are defined on "parent" mesh,
 *                  false if values are defined on post-processing mesh
 *   var_type   <-- variable's data type
 *   nt_cur_abs <-- current time step number
 *   t_cur_abs  <-- current physical time
 *   vtx_vals   <-- vertex values
 *----------------------------------------------------------------------------*/

void
cs_post_write_vertex_var(int              mesh_id,
                         const char      *var_name,
                         cs_int_t         var_dim,
                         bool             interlace,
                         bool             use_parent,
                         cs_post_type_t   var_type,
                         cs_int_t         nt_cur_abs,
                         cs_real_t        t_cur_abs,
                         const void      *vtx_vals);

/*----------------------------------------------------------------------------
 * Update references to parent mesh of post-processing meshes in case of
 * computational mesh cell renumbering.
 *
 * This function may be called only once, after possible renumbering of cells,
 * to update existing post-processing meshes. Post-processing meshes defined
 * after renumbering will automatically be based upon the new numbering,
 * so this function will not need to be called again.
 *
 * parameters:
 *   init_cell_num <-- initial cell numbering (1 to n, new -> old)
 *----------------------------------------------------------------------------*/

void
cs_post_renum_cells(const cs_int_t  init_cell_num[]);

/*----------------------------------------------------------------------------
 * Update references to parent mesh of post-processing meshes in case of
 * computational mesh interior and/or boundary faces renumbering.
 *
 * This function may be called only once, after possible renumbering of faces,
 * to update existing post-processing meshes. Post-processing meshes defined
 * after renumbering will automatically be based upon the new numbering,
 * so this function will not need to be called again.
 *
 * parameters:
 *   init_i_face_num <-- initial interior numbering (1 to n, new -> old)
 *   init_b_face_num <-- initial boundary numbering (1 to n, new -> old)
 *----------------------------------------------------------------------------*/

void
cs_post_renum_faces(const cs_int_t  init_i_face_num[],
                    const cs_int_t  init_b_face_num[]);

/*----------------------------------------------------------------------------
 * Initialize post-processing writers
 *----------------------------------------------------------------------------*/

void
cs_post_init_writers(void);

/*----------------------------------------------------------------------------
 * Initialize main post-processing meshes
 *
 * The check_flag variable is a mask, used for additionnal post-processing:
 *
 *  - If (check_flag & 1), volume submeshes are output by groups if more
 *    than one group is present and the default writer uses the EnSight format.
 *
 *  - If (check_flag & 2), boundary submeshes are output by groups if more
 *    than one group is present and the default writer uses the EnSight format.
 *
 * Note that all alias-type post-processing meshes and the meshes they
 * relate to should have been defined before calling this function, so it is
 * recommended that user-defined post-processing meshes be defined before
 * calling this function, though specific "automatic" meshes (for example
 * those related to couplings) may be defined between this call and a
 * time loop.
 *
 * parameters:
 *   check_flag <-- mask used for additional output
 *----------------------------------------------------------------------------*/

void
cs_post_init_meshes(int check_mask);

/*----------------------------------------------------------------------------
 * Destroy all structures associated with post-processing
 *----------------------------------------------------------------------------*/

void
cs_post_finalize(void);

/*----------------------------------------------------------------------------
 * Postprocess free (isolated) faces of the current global mesh
 *----------------------------------------------------------------------------*/

void
cs_post_add_free_faces(void);

/*----------------------------------------------------------------------------
 * Initialize post-processing writer with same format and associated
 * options as default writer, but no time dependency, intended to
 * troubleshoot errors.
 *----------------------------------------------------------------------------*/

void
cs_post_init_error_writer(void);

/*----------------------------------------------------------------------------
 * Initialize post-processing writer with same format and associated
 * options as default writer, but no time dependency, and associate
 * and output global volume mesh.
 *
 * This is intended to help troubleshoot errors using fields based
 * on cells.
 *
 * returns:
 *   id of error output mesh (< 0), or 0 if all writers are deactivated
 *----------------------------------------------------------------------------*/

int
cs_post_init_error_writer_cells(void);

/*----------------------------------------------------------------------------
 * Register a processing of a time-dependent variable to the call to PSTVAR.
 *
 * The instance identifier associated with the function allows registering
 * the same function several times, with a diferent identifier allowing the
 * function to select a specific operation or data.
 *
 * parameters:
 *   function    <-- function to register
 *   instance_id <-- instance id associated with this registration
 *----------------------------------------------------------------------------*/

void
cs_post_add_time_dep_var(cs_post_time_dep_var_t  *function,
                         int                      instance_id);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_POST_H__ */
