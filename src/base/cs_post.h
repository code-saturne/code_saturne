#ifndef __CS_POST_H__
#define __CS_POST_H__

/*============================================================================
 * Post-processing management
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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

#include "cs_base.h"
#include "cs_probe.h"
#include "cs_time_step.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*
 * Output type masks
 */

#define CS_POST_ON_LOCATION          (1 << 0)  /* postprocess variables
                                                  on their base location
                                                  (volume for variables) */
#define CS_POST_BOUNDARY_NR          (1 << 1)  /* postprocess boundary
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

/*----------------------------------------------------------------------------
 * Function pointer to elements selection definition.
 *
 * Each function of this sort may be used to select a given type of element,
 * usually cells, interior faces, boundary faces, or particles.
 *
 * If non-empty and not containing all elements, a list of elements of the
 * main mesh should be allocated (using BFT_MALLOC) and defined by this
 * function when called. This list's lifecycle is then managed by the
 * postprocessing subsystem.
 *
 * Note: if the input pointer is non-NULL, it must point to valid data
 * when the selection function is called, so either:
 * - that value or structure should not be temporary (i.e. local);
 * - post-processing output must be ensured using cs_post_write_meshes()
 *   with a fixed-mesh writer before the data pointed to goes out of scope;
 *
 * parameters:
 *   input    <-> pointer to optional (untyped) value or structure.
 *   n_elts   --> number of selected elements.
 *   elt_list --> list of selected elements (0 to n-1 numbering).
 *----------------------------------------------------------------------------*/

typedef void
(cs_post_elt_select_t) (void        *input,
                        cs_lnum_t   *n_elts,
                        cs_lnum_t  **elt_list);

/*----------------------------------------------------------------------------
 * Function pointer associated with a specific post-processing output.
 *
 * Such functions are registered using the cs_post_add_time_dep_vars(),
 * and all registered functions are automatically called by
 * cs_post_write_vars().
 *
 * Note: if the input pointer is non-NULL, it must point to valid data
 * when the output function is called, so either:
 * - that value or structure should not be temporary (i.e. local);
 * - post-processing output must be ensured using cs_post_write_var()
 *   or similar before the data pointed to goes out of scope.
 *
 * parameters:
 *   input <-> pointer to optional (untyped) value or structure.
 *   ts    <-- time step status structure, or NULL
 *----------------------------------------------------------------------------*/

typedef void
(cs_post_time_dep_output_t) (void                  *input,
                             const cs_time_step_t  *ts);

/*----------------------------------------------------------------------------
 * Function pointer associated with a specific post-processing output
 * on multiple meshes.
 *
 * Such functions are registered using the cs_post_add_time_mesh_dep_vars(),
 * and all registered functions are automatically called by
 * cs_post_write_vars().
 *
 * Note: if the input pointer is non-NULL, it must point to valid data
 * when the output function is called, so either:
 * - that value or structure should not be temporary (i.e. local);
 * - post-processing output must be ensured using cs_post_write_var()
 *   or similar before the data pointed to goes out of scope.
 *
 * parameters:
 *   input       <-> pointer to optional (untyped) value or structure.
 *   mesh_id     <-- id of the output mesh for the current call
 *   cat_id      <-- category id of the output mesh for the current call
 *   ent_flag    <-- indicate global presence of cells (ent_flag[0]), interior
 *                   faces (ent_flag[1]), boundary faces (ent_flag[2]),
 *                   particles (ent_flag[3]) or probes (ent_flag[4])
 *   n_cells     <-- local number of cells of post_mesh
 *   n_i_faces   <-- local number of interior faces of post_mesh
 *   n_b_faces   <-- local number of boundary faces of post_mesh
 *   cell_list   <-- list of cells (1 to n) of post-processing mesh
 *   i_face_list <-- list of interior faces (1 to n) of post-processing mesh
 *   b_face_list <-- list of boundary faces (1 to n) of post-processing mesh
 *   ts          <-- time step status structure, or NULL
 *----------------------------------------------------------------------------*/

typedef void
(cs_post_time_mesh_dep_output_t) (void                  *input,
                                  int                    mesh_id,
                                  int                    cat_id,
                                  int                    ent_flag[5],
                                  cs_lnum_t              n_cells,
                                  cs_lnum_t              n_i_faces,
                                  cs_lnum_t              n_b_faces,
                                  const cs_lnum_t        cell_list[],
                                  const cs_lnum_t        i_face_list[],
                                  const cs_lnum_t        b_face_list[],
                                  const cs_time_step_t  *ts);

/*=============================================================================
 * Global variables
 *============================================================================*/

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
 * Configure the post-processing output so that mesh connectivity
 * may be automatically updated.
 *
 * This is done for meshes defined using selection criteria or functions.
 * The behavior of Lagrangian meshes is unchanged.
 *
 * To be effective, this function should be called before defining
 * postprocessing meshes.
 *----------------------------------------------------------------------------*/

void
cs_post_set_changing_connectivity(void);

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
                      int                     frequency_n,
                      double                  frequency_t);

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
 * Define a volume post-processing mesh using a selection function.
 *
 * The selection may be updated over time steps if both the time_varying
 * flag is set to true and the mesh is only associated with writers defined
 * with the FVM_WRITER_TRANSIENT_CONNECT option.
 *
 * Note: if the cell_select_input pointer is non-NULL, it must point
 * to valid data when the selection function is called, so either:
 * - that value or structure should not be temporary (i.e. local);
 * - post-processing output must be ensured using cs_post_write_meshes()
 *   with a fixed-mesh writer before the data pointed to goes out of scope;
 *
 * parameters:
 *   mesh_id           <-- id of mesh to define (< 0 reserved, > 0 for user)
 *   mesh_name         <-- associated mesh name
 *   cell_select_func  <-- pointer to cells selection function
 *   cell_select_input <-> pointer to optional input data for the cell
 *                         selection function, or NULL
 *   time_varying      <-- if true, try to redefine mesh at each output time
 *   add_groups        <-- if true, add group information if present
 *   auto_variables    <-- if true, automatic output of main variables
 *   n_writers         <-- number of associated writers
 *   writer_ids        <-- ids of associated writers
 *----------------------------------------------------------------------------*/

void
cs_post_define_volume_mesh_by_func(int                    mesh_id,
                                   const char            *mesh_name,
                                   cs_post_elt_select_t  *cell_select_func,
                                   void                  *cell_select_input,
                                   bool                   time_varying,
                                   bool                   add_groups,
                                   bool                   auto_variables,
                                   int                    n_writers,
                                   const int              writer_ids[]);

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
 * Define a surface post-processing mesh using selection functions.
 *
 * The selection may be updated over time steps if both the time_varying
 * flag is set to true and the mesh is only associated with writers defined
 * with the FVM_WRITER_TRANSIENT_CONNECT option.
 *
 * Note: if i_face_select_input or b_face_select_input pointer is non-NULL,
 * it must point to valid data when the selection function is called,
 * so either:
 * - that value or structure should not be temporary (i.e. local);
 * - post-processing output must be ensured using cs_post_write_meshes()
 *   with a fixed-mesh writer before the data pointed to goes out of scope;
 *
 * parameters:
 *   mesh_id             <-- id of mesh to define (< 0 reserved, > 0 for user)
 *   mesh_name           <-- associated mesh name
 *   i_face_select_func  <-- pointer to interior faces selection function
 *   b_face_select_func  <-- pointer to boundary faces selection function
 *   i_face_select_input <-> pointer to optional input data for the interior
 *                           faces selection function, or NULL
 *   b_face_select_input <-> pointer to optional input data for the boundary
 *                           faces selection function, or NULL
 *   time_varying        <-- if true, try to redefine mesh at each output time
 *   add_groups          <-- if true, add group information if present
 *   auto_variables      <-- if true, automatic output of main variables
 *   n_writers           <-- number of associated writers
 *   writer_ids          <-- ids of associated writers
 *----------------------------------------------------------------------------*/

void
cs_post_define_surface_mesh_by_func(int                    mesh_id,
                                    const char            *mesh_name,
                                    cs_post_elt_select_t  *i_face_select_func,
                                    cs_post_elt_select_t  *b_face_select_func,
                                    void                  *i_face_select_input,
                                    void                  *b_face_select_input,
                                    bool                   time_varying,
                                    bool                   add_groups,
                                    bool                   auto_variables,
                                    int                    n_writers,
                                    const int              writer_ids[]);

/*----------------------------------------------------------------------------
 * Define a particles post-processing mesh.
 *
 * Such a mesh is always time-varying, and will only be output by writers
 * defined with the FVM_WRITER_TRANSIENT_CONNECT option.
 *
 * If the trajectory_mode argument is set to true, this logic is reversed,
 * and output will only occur for writers defined with the
 * FVM_WRITER_FIXED_MESH option. In this case, a submesh consisting of
 * trajectory segments for the current time step will be added to
 * the output at each output time step.
 *
 * parameters:
 *   mesh_id        <-- id of mesh to define (< 0 reserved, > 0 for user)
 *   mesh_name      <-- associated mesh name
 *   cell_criteria  <-- selection criteria for cells containing particles,
 *                      or NULL.
 *   density        <-- fraction of the particles in the selected area
 *                      which should be output (0 < density <= 1)
 *   trajectory     <-- if true, activate trajectory mode
 *   auto_variables <-- if true, automatic output of main variables
 *   n_writers      <-- number of associated writers
 *   writer_ids     <-- ids of associated writers
 *----------------------------------------------------------------------------*/

void
cs_post_define_particles_mesh(int          mesh_id,
                              const char  *mesh_name,
                              const char  *cell_criteria,
                              double       density,
                              bool         trajectory,
                              bool         auto_variables,
                              int          n_writers,
                              const int    writer_ids[]);

/*----------------------------------------------------------------------------
 * Define a particles post-processing mesh using a selection function.
 *
 * The selection may be updated over time steps.
 *
 * Such a mesh is always time-varying, and will only be output by writers
 * defined with the FVM_WRITER_TRANSIENT_CONNECT option.
 *
 * If the trajectory_mode argument is set to true, this logic is reversed,
 * and output will only occur for writers defined with the
 * FVM_WRITER_FIXED_MESH option. In this case, a submesh consisting of
 * trajectory segments for the current time step will be added to
 * the output at each output time step.
 *
 * Note: if the p_select_input pointer is non-NULL, it must point
 * to valid data when the selection function is called, so
 * that value or structure should not be temporary (i.e. local);
 *
 * parameters:
 *   mesh_id        <-- id of mesh to define (< 0 reserved, > 0 for user)
 *   mesh_name      <-- associated mesh name
 *   p_select_func  <-- pointer to particles selection function
 *   p_select_input <-> pointer to optional input data for the particles
 *                      selection function, or NULL
 *   density        <-- fraction of the particles in the selected area
 *                      which should be output (0 < density <= 1)
 *   trajectory     <-- if true, activate trajectory mode
 *   auto_variables <-- if true, automatic output of main variables
 *   n_writers      <-- number of associated writers
 *   writer_ids     <-- ids of associated writers
 *----------------------------------------------------------------------------*/

void
cs_post_define_particles_mesh_by_func(int                    mesh_id,
                                      const char            *mesh_name,
                                      cs_post_elt_select_t  *p_select_func,
                                      void                  *p_select_input,
                                      bool                   trajectory,
                                      bool                   auto_variables,
                                      int                    n_writers,
                                      const int              writer_ids[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a post-processing mesh for probes (set of probes should have
 *        been already defined)
 *
 * \param[in] mesh_id        id of mesh to define (<0 reserved, >0 for user)
 * \param[in] pset           pointer to a cs_probe_set_t structure
 * \param[in] time_varying   true if probe coords may change during computation
 * \param[in] is_profile     true if probe set is related to a profile
 * \param[in] on_boundary    true if probes are located on boundary
 * \param[in] auto_variable  true if the set of variables to post is predefined
 * \param[in] n_writers      number of associated writers
 * \param[in] writer_ids     ids of associated writers
 */
/*----------------------------------------------------------------------------*/

void
cs_post_define_probe_mesh(int                    mesh_id,
                          cs_probe_set_t        *pset,
                          bool                   time_varying,
                          bool                   is_profile,
                          bool                   on_boundary,
                          bool                   auto_variable,
                          int                    n_writers,
                          const int              writer_ids[]);

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
 *   mesh_id      <-- id of edges mesh to create (< 0 reserved, > 0 for user)
 *   base_mesh_id <-- id of existing mesh (< 0 reserved, > 0 for user)
 *   n_writers    <-- number of associated writers
 *   writer_ids   <-- ids of associated writers
 *----------------------------------------------------------------------------*/

void
cs_post_define_edges_mesh(int        mesh_id,
                          int        base_mesh_id,
                          int        n_writers,
                          const int  writer_ids[]);

/*----------------------------------------------------------------------------
 * Get a postprocessing meshes entity presence flag.
 *
 * This flag is an array of 3 integers, indicating the presence of elements
 * of given types on at least one subdomain (i.e. rank):
 *   0: presence of cells
 *   1: presence of interior faces
 *   2: presence of boundary faces
 *
 * parameters:
 *   mesh_id <-- postprocessing mesh id
 *
 * returns:
 *   pointer to entity presence flag
 *----------------------------------------------------------------------------*/

const int *
cs_post_mesh_get_ent_flag(int  mesh_id);

/*----------------------------------------------------------------------------
 * Get a postprocessing mesh's number of cells.
 *
 * parameters:
 *   mesh_id <-- postprocessing mesh id
 *
 * returns:
 *   number of cells of postprocessing mesh.
 *----------------------------------------------------------------------------*/

cs_lnum_t
cs_post_mesh_get_n_cells(int  mesh_id);

/*----------------------------------------------------------------------------
 * Get a postprocessing mesh's list of cells.
 *
 * The array of cell ids must be of at least size
 * cs_post_mesh_get_n_cells(mesh_id).
 *
 * parameters:
 *   mesh_id  <-- postprocessing mesh id
 *   cell_ids --> array of associated cell ids (0 to n-1 numbering,
 *                relative to main mesh)
 *----------------------------------------------------------------------------*/

void
cs_post_mesh_get_cell_ids(int         mesh_id,
                          cs_lnum_t  *cell_ids);

/*----------------------------------------------------------------------------
 * Get a postprocessing mesh's number of interior faces.
 *
 * parameters:
 *   mesh_id <-- postprocessing mesh id
 *
 * returns:
 *   number of cells of postprocessing mesh.
 *----------------------------------------------------------------------------*/

cs_lnum_t
cs_post_mesh_get_n_i_faces(int  mesh_id);

/*----------------------------------------------------------------------------
 * Get a postprocessing mesh's list of boundary faces.
 *
 * The array of boundary face ids must be of at least size
 * cs_post_mesh_get_n_b_faces(mesh_id).
 *
 * parameters:
 *   mesh_id    <-- postprocessing mesh id
 *   i_face_ids --> array of associated interior faces ids
 *                  (0 to n-1 numbering, relative to main mesh)
 *----------------------------------------------------------------------------*/

void
cs_post_mesh_get_i_face_ids(int        mesh_id,
                            cs_lnum_t  i_face_ids[]);

/*----------------------------------------------------------------------------
 * Get a postprocessing mesh's number of boundary faces
 *
 * parameters:
 *   mesh_id <-- postprocessing mesh id
 *
 * returns:
 *   number of cells of postprocessing mesh.
 *----------------------------------------------------------------------------*/

cs_lnum_t
cs_post_mesh_get_n_b_faces(int  mesh_id);

/*----------------------------------------------------------------------------
 * Get a postprocessing mesh's list of boundary faces.
 *
 * The array of boundary face ids must be of at least size
 * cs_post_mesh_get_n_b_faces(mesh_id).
 *
 * parameters:
 *   mesh_id    <-- postprocessing mesh id
 *   b_face_ids --> array of associated boundary faces ids
 *                  (0 to n-1 numbering, relative to main mesh)
 *----------------------------------------------------------------------------*/

void
cs_post_mesh_get_b_face_ids(int        mesh_id,
                            cs_lnum_t  b_face_ids[]);

/*----------------------------------------------------------------------------
 * Set whether postprocessing mesh's parallel domain should be output.
 *
 * parameters:
 *   mesh_id     <-- id of mesh to remove
 *   post_domain <- true if parallel domain should be output,
 *                  false otherwise.
 *----------------------------------------------------------------------------*/

void
cs_post_mesh_set_post_domain(int   mesh_id,
                             bool  post_domain);

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
 * Return a pointer to the FVM writer associated to a writer_id.
 *
 * parameters:
 *   writer_id <-- associated writer id
 *
 * returns:
 *  a pointer to a fvm_writer_t structure
 *----------------------------------------------------------------------------*/

fvm_writer_t *
cs_post_get_writer(int  writer_id);

/*----------------------------------------------------------------------------
 * Return time dependency associated to a writer_id.
 *
 * parameters:
 *   writer_id <-- associated writer id
 *
 * returns:
 *   associated writer's time dependency
 *----------------------------------------------------------------------------*/

fvm_writer_time_dep_t
cs_post_get_writer_time_dep(int  writer_id);

/*----------------------------------------------------------------------------
 * Add an activation time step for a specific writer or for all writers.
 *
 * If a negative value is provided, a previously added activation time
 * step matching that absolute value will be removed, if present.
 *
 * parameters:
 *   writer_id <-- writer id, or 0 for all writers
 *   nt        <-- time step value to add (or remove)
 *----------------------------------------------------------------------------*/

void
cs_post_add_writer_t_step(int  writer_id,
                          int  nt);

/*----------------------------------------------------------------------------
 * Add an activation time value for a specific writer or for all writers.
 *
 * If a negative value is provided, a previously added activation time
 * step matching that absolute value will be removed, if present.
 *
 * parameters:
 *   writer_id <-- writer id, or 0 for all writers
 *   t         <-- time value to add (or remove)
 *----------------------------------------------------------------------------*/

void
cs_post_add_writer_t_value(int     writer_id,
                           double  t);

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
cs_post_modify_mesh(int        mesh_id,
                    cs_lnum_t  n_cells,
                    cs_lnum_t  n_i_faces,
                    cs_lnum_t  n_b_faces,
                    cs_lnum_t  cell_list[],
                    cs_lnum_t  i_face_list[],
                    cs_lnum_t  b_face_list[]);

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
 * Update "active" or "inactive" flag of writers based on the time step.
 *
 * Writers are activated if their output frequency is a divisor of the
 * current time step, or if their optional time step and value output lists
 * contain matches for the current time step.
 *
 * parameters:
 *   ts <-- time step status structure
 *----------------------------------------------------------------------------*/

void
cs_post_activate_by_time_step(const cs_time_step_t  *ts);

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
 *   ts <-- time step status structure
 *----------------------------------------------------------------------------*/

void
cs_post_write_meshes(const cs_time_step_t  *ts);

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
 *   cel_vals    <-- cell values
 *   i_face_vals <-- interior face values
 *   b_face_vals <-- boundary face values
 *   ts          <-- time step status structure, or NULL
 *----------------------------------------------------------------------------*/

void
cs_post_write_var(int                    mesh_id,
                  const char            *var_name,
                  int                    var_dim,
                  bool                   interlace,
                  bool                   use_parent,
                  cs_post_type_t         var_type,
                  const void            *cel_vals,
                  const void            *i_face_vals,
                  const void            *b_face_vals,
                  const cs_time_step_t  *ts);

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
 *   vtx_vals   <-- vertex values
 *   ts          <-- time step status structure, or NULL
 *----------------------------------------------------------------------------*/

void
cs_post_write_vertex_var(int                    mesh_id,
                         const char            *var_name,
                         int                    var_dim,
                         bool                   interlace,
                         bool                   use_parent,
                         cs_post_type_t         var_type,
                         const void            *vtx_vals,
                         const cs_time_step_t  *ts);

/*----------------------------------------------------------------------------
 * Output an existing lagrangian particle attribute at particle
 * positions or trajectory endpoints of a particle mesh using
 * associated writers.
 *
 * parameters:
 *   mesh_id      <-- id of associated mesh
 *   attr         <-- associated particle attribute id
 *   var_name     <-- name of variable to output
 *   component_id <-- if -1 : extract the whole attribute
 *                    if >0 : id of the component to extract
 *   ts           <-- time step status structure, or NULL
 *----------------------------------------------------------------------------*/

void
cs_post_write_particle_values(int                    mesh_id,
                              int                    attr_id,
                              const char            *var_name,
                              int                    component_id,
                              const cs_time_step_t  *ts);

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
 *   init_cell_num <-- initial cell numbering (new -> old)
 *----------------------------------------------------------------------------*/

void
cs_post_renum_cells(const cs_lnum_t  init_cell_num[]);

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
 *   init_i_face_num <-- initial interior numbering (new -> old)
 *   init_b_face_num <-- initial boundary numbering (new -> old)
 *----------------------------------------------------------------------------*/

void
cs_post_renum_faces(const cs_lnum_t  init_i_face_num[],
                    const cs_lnum_t  init_b_face_num[]);

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
 * Loop on post-processing meshes to output variables.
 *
 * This handles all default fields output, as well as all
 * registred output functions.
 *
 * parameters:
 *   ts <-- time step status structure, or NULL
 *----------------------------------------------------------------------------*/

void
cs_post_write_vars(const cs_time_step_t  *ts);

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
 * Register a processing of time-dependent variables to the call to
 * cs_post_write_vars().
 *
 * Note: if the input pointer is non-NULL, it must point to valid data
 * when the output function is called, so either:
 * - that value or structure should not be temporary (i.e. local);
 * - post-processing output must be ensured using cs_post_write_var()
 *   or similar before the data pointed to goes out of scope.
 *
 * parameters:
 *   function <-- function to register
 *   input    <-> pointer to optional (untyped) value or structure.
 *----------------------------------------------------------------------------*/

void
cs_post_add_time_dep_output(cs_post_time_dep_output_t  *function,
                            void                       *input);

/*----------------------------------------------------------------------------
 * Register a processing of time-dependent variables than can be output
 * on different meshes to the call to cs_post_write_vars().
 *
 * Note: if the input pointer is non-NULL, it must point to valid data
 * when the output function is called, so either:
 * - that value or structure should not be temporary (i.e. local);
 * - post-processing output must be ensured using cs_post_write_var()
 *   or similar before the data pointed to goes out of scope.
 *
 * parameters:
 *   function <-- function to register
 *   input    <-> pointer to optional (untyped) value or structure.
 *----------------------------------------------------------------------------*/

void
cs_post_add_time_mesh_dep_output(cs_post_time_mesh_dep_output_t  *function,
                                 void                            *input);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_POST_H__ */
