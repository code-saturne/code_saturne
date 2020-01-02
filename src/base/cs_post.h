#ifndef __CS_POST_H__
#define __CS_POST_H__

/*============================================================================
 * Post-processing management
 *============================================================================*/

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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_nodal.h"
#include "fvm_writer.h"

#include "cs_base.h"
#include "cs_interpolate.h"
#include "cs_probe.h"
#include "cs_time_step.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/* Output type masks */

#define CS_POST_ON_LOCATION          (1 << 0)  /* postprocess variables
                                                  on their base location
                                                  (volume for variables) */
#define CS_POST_BOUNDARY_NR          (1 << 1)  /* postprocess boundary
                                                  without reconstruction */

#define CS_POST_MONITOR              (1 << 2)  /* monitor variables */

/* Default writer ids and filters */

#define CS_POST_WRITER_ALL_ASSOCIATED  0       /* all associated writers */

#define CS_POST_WRITER_DEFAULT        -1       /* default visualisation */
#define CS_POST_WRITER_ERRORS         -2       /* error visualisation */
#define CS_POST_WRITER_PARTICLES      -3       /* particle visualisation */
#define CS_POST_WRITER_TRAJECTORIES   -4       /* trajectories visualisation */
#define CS_POST_WRITER_PROBES         -5       /* probe monitoring */
#define CS_POST_WRITER_PROFILES       -6       /* profiles */
#define CS_POST_WRITER_HISTOGRAMS     -7       /* histograms */

/* Default mesh ids */

#define CS_POST_MESH_VOLUME           -1       /* volume mesh output */
#define CS_POST_MESH_BOUNDARY         -2       /* boundary mesh output */
#define CS_POST_MESH_PARTICLES        -3       /* particle output */
#define CS_POST_MESH_TRAJECTORIES     -4       /* particle output */
#define CS_POST_MESH_PROBES           -5       /* probes output */

/* Additional categories (no associated default mesh) */

#define CS_POST_MESH_SURFACE         -12       /* surface (boundary and/or
                                                  interior) mesh */

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
 *   cell_ids    <-- list of cells (0 to n-1) of post-processing mesh
 *   i_face_ids  <-- list of interior faces (0 to n-1) of post-processing mesh
 *   b_face_ids  <-- list of boundary faces (0 to n-1) of post-processing mesh
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
                                  const cs_lnum_t        cell_ids[],
                                  const cs_lnum_t        i_face_ids[],
                                  const cs_lnum_t        b_face_ids[],
                                  const cs_time_step_t  *ts);

/*=============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a writer; this objects manages a case's name, directory,
 *        and format, as well as associated mesh's time dependency, and the
 *        default output frequency for associated variables.
 *
 * This function must be called before the time loop. If a writer with a
 * given id is defined multiple times, the last definition supercedes the
 * previous ones.
 *
 * Current reserved ids are the following: CS_POST_WRITER_DEFAULT
 * for main/default output, CS_POST_WRITER_ERRORS for error visualization,
 * CS_POST_WRITER_PROBES for main probes, CS_POST_WRITER_PARTICLES for
 * particles, CS_POST_WRITER_TRAJECTORIES for trajectories. Other negative
 * ids may be dynamically reserved by the code depending on options.
 * Positive ids identify user-defined writers.
 *
 * \warning depending on the chosen format, the \em case_name may be
 * shortened (maximum number of characters: 32 for \em MED, 19 for \em EnSight,
 * or modified automatically (white-space or forbidden characters will be
 * replaced by "_").
 *
 * The \c \b format_name argument is used to choose the output format, and the
 * following values are allowed (assuming the matching
 * support was built):
 *
 * - \c \b EnSight \c \b Gold (\c \b EnSight also accepted)
 * - \c \b MED
 * - \c \b CGNS
 * - \c \b CCM (only for the full volume and boundary meshes)
 * - \c \b Catalyst (in-situ visualization)
 * - \c \b MEDCoupling (in-memory structure, to be used from other code)
 * - \c \b plot (comma or whitespace separated 2d plot files)
 * - \c \b time_plot (comma or whitespace separated time plot files)
 *
 * The format name is case-sensitive, so \c \b ensight or \c \b cgns are also valid.
 *
 * The optional \c \b fmt_opts character string contains a list of options related
 * to the format, separated by spaces or commas; these options include:
 *
 * - \c \b binary for a binary format version (default)
 * - \c \b big_endian to force outputs to be in \c \b big-endian mode
 *         (for \c \b EnSight).
 * - \c \b text for a text format version (for \c \b EnSight).
 * - \c \b adf for ADF file type (for \c \b CGNS).
 * - \c \b hdf5 for HDF5 file type (for \c \b CGNS, normally the default if
 *         HDF5 support is available).
 * - \c \b discard_polygons to prevent from exporting faces with more than
 *         four edges (which may not be recognized by some post-processing
 *         tools); such faces will therefore not appear in the post-processing
 *         mesh.
 * - \c \b discard_polyhedra to prevent from exporting elements which are
 *         neither tetrahedra, prisms, pyramids nor hexahedra (which may not
 *         be recognized by some post-processing tools); such elements will
 *         therefore not appear in the post-processing mesh.
 * - \c \b divide_polygons to divide faces with more than four edges into
 *         triangles, so that any post-processing tool can recognize them.
 * - \c \b divide_polyhedra} to divide elements which are neither tetrahedra,
 *         prisms, pyramids nor hexahedra into simpler elements (tetrahedra and
 *         pyramids), so that any post-processing tool can recognize them.
 * - \c \b separate_meshes to multiple meshes and associated fields to
 *         separate outputs.
 *
 * Note that the white-spaces in the beginning or in the end of the
 * character strings given as arguments here are suppressed automatically.
 *
 * \param[in]  writer_id        id of writer to create. (< 0 reserved,
 *                              > 0 for user); eveb for reserved ids,
 *                              the matching writer's options
 *                              may be redifined by calls to this function
 * \param[in]  case_name        associated case name
 * \param[in]  dir_name         associated directory name
 * \param[in]  fmt_name         associated format name
 * \param[in]  fmt_opts         associated format options string
 * \param[in]  time_dep         \ref FVM_WRITER_FIXED_MESH if mesh definitions
 *                              are fixed, \ref FVM_WRITER_TRANSIENT_COORDS if
 *                              coordinates change,
 *                              \ref FVM_WRITER_TRANSIENT_CONNECT if
 *                              connectivity changes
 * \param[in]  output_at_start  force output at calculation start if true
 * \param[in]  output_at_end    force output at calculation end if true
 * \param[in]  frequency_n      default output frequency in time-steps, or < 0
 * \param[in]  frequency_t      default output frequency in seconds, or < 0
 *                              (has priority over frequency_n)
 */
/*----------------------------------------------------------------------------*/

void
cs_post_define_writer(int                     writer_id,
                      const char             *case_name,
                      const char             *dir_name,
                      const char             *fmt_name,
                      const char             *fmt_opts,
                      fvm_writer_time_dep_t   time_dep,
                      bool                    output_at_start,
                      bool                    output_at_end,
                      int                     frequency_n,
                      double                  frequency_t);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a volume post-processing mesh.
 *
 * \param[in]  mesh_id         id of mesh to define
 *                             (< 0 reserved, > 0 for user)
 * \param[in]  mesh_name       associated mesh name
 * \param[in]  cell_criteria   selection criteria for cells
 * \param[in]  add_groups      if true, add group information if present
 * \param[in]  auto_variables  if true, automatic output of main variables
 * \param[in]  n_writers       number of associated writers
 * \param[in]  writer_ids      ids of associated writers
 */
/*----------------------------------------------------------------------------*/

void
cs_post_define_volume_mesh(int          mesh_id,
                           const char  *mesh_name,
                           const char  *cell_criteria,
                           bool         add_groups,
                           bool         auto_variables,
                           int          n_writers,
                           const int    writer_ids[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a volume post-processing mesh using a selection function.
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
 * \param[in]  mesh_id            id of mesh to define
 *                                (< 0 reserved, > 0 for user)
 * \param[in]  mesh_name          associated mesh name
 * \param[in]  cell_select_func   pointer to cells selection function
 * \param[in]  cell_select_input  pointer to optional input data for the cell
 *                                selection function, or NULL
 * \param[in]  time_varying       if true, try to redefine mesh at each
 *                                output time
 * \param[in]  add_groups         if true, add group information if present
 * \param[in]  auto_variables     if true, automatic output of main variables
 * \param[in]  n_writers          number of associated writers
 * \param[in]  writer_ids         ids of associated writers
 */
/*----------------------------------------------------------------------------*/

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

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a surface post-processing mesh.
 *
 * \param[in]  mesh_id          id of mesh to define
 *                              (< 0 reserved, > 0 for user)
 * \param[in]  mesh_name        associated mesh name
 * \param[in]  i_face_criteria  selection criteria for interior faces
 * \param[in]  b_face_criteria  selection criteria for boundary faces
 * \param[in]  add_groups       if true, add group information if present
 * \param[in]  auto_variables   if true, automatic output of main variables
 * \param[in]  n_writers        number of associated writers
 * \param[in]  writer_ids       ids of associated writers
 */
/*----------------------------------------------------------------------------*/

void
cs_post_define_surface_mesh(int          mesh_id,
                            const char  *mesh_name,
                            const char  *i_face_criteria,
                            const char  *b_face_criteria,
                            bool         add_groups,
                            bool         auto_variables,
                            int          n_writers,
                            const int    writer_ids[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a surface post-processing mesh using selection functions.
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
 * \param[in]  mesh_id              id of mesh to define
 *                                  (< 0 reserved, > 0 for user)
 * \param[in]  mesh_name            associated mesh name
 * \param[in]  i_face_select_func   pointer to interior faces selection function
 * \param[in]  b_face_select_func   pointer to boundary faces selection function
 * \param[in]  i_face_select_input  pointer to optional input data for the
 *                                  interior faces selection function, or NULL
 * \param[in]  b_face_select_input  pointer to optional input data for the
 *                                  boundary faces selection function, or NULL
 * \param[in]  time_varying         if true, try to redefine mesh at each
 *                                  output time
 * \param[in]  add_groups           if true, add group information if present
 * \param[in]  auto_variables       if true, automatic output of main variables
 * \param[in]  n_writers            number of associated writers
 * \param[in]  writer_ids          ids of associated writers
 */
/*----------------------------------------------------------------------------*/

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

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a particles post-processing mesh.
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
 * \param[in]  mesh_id         id of mesh to define
 *                             (< 0 reserved, > 0 for user)
 * \param[in]  mesh_name       associated mesh name
 * \param[in]  cell_criteria   selection criteria for cells containing
 *                             particles, or NULL.
 * \param[in]  density         fraction of the particles in the selected area
 *                             which should be output (0 < density <= 1)
 * \param[in]  trajectory      if true, activate trajectory mode
 * \param[in]  auto_variables  if true, automatic output of main variables
 * \param[in]  n_writers       number of associated writers
 * \param[in]  writer_ids      ids of associated writers
 */
/*----------------------------------------------------------------------------*/

void
cs_post_define_particles_mesh(int          mesh_id,
                              const char  *mesh_name,
                              const char  *cell_criteria,
                              double       density,
                              bool         trajectory,
                              bool         auto_variables,
                              int          n_writers,
                              const int    writer_ids[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a particles post-processing mesh using a selection function.
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
 * \param[in]  mesh_id         id of mesh to define
 *                             (< 0 reserved, > 0 for user)
 * \param[in]  mesh_name       associated mesh name
 * \param[in]  p_select_func   pointer to particles selection function
 * \param[in]  p_select_input  pointer to optional input data for the particles
 *                             selection function, or NULL
 * \param[in]  trajectory      if true, activate trajectory mode
 * \param[in]  auto_variables  if true, automatic output of main variables
 * \param[in]  n_writers       number of associated writers
 * \param[in]  writer_ids      ids of associated writers
 */
/*----------------------------------------------------------------------------*/

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
 * \brief Create a post-processing mesh associated with an existing exportable
 * mesh representation.
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
 * \param[in]  mesh_id         id of mesh to define
 *                             (< 0 reserved, > 0 for user)
 * \param[in]  exp_mesh        mesh in exportable representation
 *                             (i.e. fvm_nodal_t)
 * \param[in]  dim_shift       nonzero if exp_mesh has been projected
 * \param[in]  transfer        if true, ownership of exp_mesh is transferred
 *                             to the post-processing mesh
 * \param[in]  auto_variables  if true, automatic output of main variables
 * \param[in]  n_writers       number of associated writers
 * \param[in]  writer_ids      ids of associated writers
 */
/*----------------------------------------------------------------------------*/

void
cs_post_define_existing_mesh(int           mesh_id,
                             fvm_nodal_t  *exp_mesh,
                             int           dim_shift,
                             bool          transfer,
                             bool          auto_variables,
                             int           n_writers,
                             const int     writer_ids[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a mesh based upon the extraction of edges from an existing mesh.
 *
 * The newly created edges have no link to their parent elements, so
 * no variable referencing parent elements may be output to this mesh,
 * whose main use is to visualize "true" face edges when polygonal faces
 * are subdivided by the writer. In this way, even highly non-convex
 * faces may be visualized correctly if their edges are overlaid on
 * the surface mesh with subdivided polygons.
 *
 * \param[in]  mesh_id       id of edges mesh to create
 *                           (< 0 reserved, > 0 for user)
 * \param[in]  base_mesh_id  id of existing mesh (< 0 reserved, > 0 for user)
 * \param[in]  n_writers     number of associated writers
 * \param[in]  writer_ids    ids of associated writers
 */
/*----------------------------------------------------------------------------*/

void
cs_post_define_edges_mesh(int        mesh_id,
                          int        base_mesh_id,
                          int        n_writers,
                          const int  writer_ids[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Associate a writer to a postprocessing mesh.
 *
 * This function must be called during the postprocessing output definition
 * stage, before any output actually occurs.
 *
 * If called with a non-existing mesh or writer id, or if the writer is
 * already associated, no setting is changed, and this function
 * returns silently.
 *
 * \param[in]  mesh_id      id of mesh to define
 *                          (< 0 reserved, > 0 for user)
 * \param[in]  writer_id    id of writer to associate
 */
/*----------------------------------------------------------------------------*/

void
cs_post_mesh_attach_writer(int  mesh_id,
                           int  writer_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief De-associate a writer from a postprocessing mesh.
 *
 * This function must be called during the postprocessing output definition
 * stage, before any output actually occurs.
 *
 * If called with a non-existing mesh or writer id, or if the writer was not
 * previously associated, no setting is changed, and this function
 * returns silently.
 *
 * \param[in]  mesh_id      id of mesh to define
 *                          (< 0 reserved, > 0 for user)
 * \param[in]  writer_id    id of writer to associate
 */
/*----------------------------------------------------------------------------*/

void
cs_post_mesh_detach_writer(int  mesh_id,
                           int  writer_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get a postprocessing meshes entity presence flag.
 *
 * This flag is an array of 5 integers, indicating the presence of elements
 * of given types on at least one subdomain (i.e. rank):
 *   0: presence of cells
 *   1: presence of interior faces
 *   2: presence of boundary faces
 *   3: presence of particles
 *   4: presence of probes
 *
 * \param[in]  mesh_id  postprocessing mesh id
 *
 * \return  pointer to entity presence flag
 */
/*----------------------------------------------------------------------------*/

const int *
cs_post_mesh_get_ent_flag(int  mesh_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get a postprocessing mesh's number of cells
 *
 * \param[in]  mesh_id  postprocessing mesh id
 *
 * \return  number of cells of postprocessing mesh.
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_post_mesh_get_n_cells(int  mesh_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get a postprocessing mesh's list of cells
 *
 * The array of cell ids must be of at least size
 * cs_post_mesh_get_n_cells(mesh_id).
 *
 * \param[in]   mesh_id   postprocessing mesh id
 * \param[out]  cell_ids  array of associated cell ids (0 to n-1 numbering,
 *                        relative to main mesh)
 */
/*----------------------------------------------------------------------------*/

void
cs_post_mesh_get_cell_ids(int         mesh_id,
                          cs_lnum_t  *cell_ids);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get a postprocessing mesh's number of interior faces
 *
 * \param[in]  mesh_id  postprocessing mesh id
 *
 * \return  number of cells of postprocessing mesh.
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_post_mesh_get_n_i_faces(int  mesh_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get a postprocessing mesh's list of boundary faces.
 *
 * The array of boundary face ids must be of at least size
 * cs_post_mesh_get_n_b_faces(mesh_id).
 *
 * \param[in]   mesh_id     postprocessing mesh id
 * \param[out]  i_face_ids  array of associated interior faces ids
 *                          (0 to n-1 numbering, relative to main mesh)
 */
/*----------------------------------------------------------------------------*/

void
cs_post_mesh_get_i_face_ids(int        mesh_id,
                            cs_lnum_t  i_face_ids[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get a postprocessing mesh's number of boundary faces
 *
 * \param[in]  mesh_id  postprocessing mesh id
 *
 * \return  number of cells of postprocessing mesh.
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_post_mesh_get_n_b_faces(int  mesh_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get a postprocessing mesh's list of boundary faces.
 *
 * The array of boundary face ids must be of at least size
 * cs_post_mesh_get_n_b_faces(mesh_id).
 *
 * \param[in]   mesh_id     postprocessing mesh id
 * \param[out]  b_face_ids  array of associated boundary faces ids
 *                          (0 to n-1 numbering, relative to main mesh)
 */
/*----------------------------------------------------------------------------*/

void
cs_post_mesh_get_b_face_ids(int        mesh_id,
                            cs_lnum_t  b_face_ids[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get a postprocessing mesh's number of vertices
 *
 * \param[in]  mesh_id  postprocessing mesh id
 *
 * \return  number of vertices of postprocessing mesh.
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_post_mesh_get_n_vertices(int  mesh_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get a postprocessing mesh's list of vertices
 *
 * The array of vertex ids must be of at least size
 * cs_post_mesh_get_n_vertices(mesh_id).
 *
 * \param[in]   mesh_id     postprocessing mesh id
 * \param[out]  vertex_ids  array of associated vertex ids (0 to n-1 numbering,
 *                          relative to main mesh)
 */
/*----------------------------------------------------------------------------*/

void
cs_post_mesh_get_vertex_ids(int         mesh_id,
                            cs_lnum_t  *vertex_ids);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set whether postprocessing mesh's parallel domain should be output.
 *
 * \param[in]  mesh_id      postprocessing mesh id
 * \param[in]  post_domain  true if parallel domain should be output,
 *                          false otherwise.
 */
/*----------------------------------------------------------------------------*/

void
cs_post_mesh_set_post_domain(int   mesh_id,
                             bool  post_domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Remove a post-processing mesh.
 *
 * No further post-processing output will be allowed on this mesh,
 * so the associated structures may be freed.
 *
 * A post-processing mesh that has been associated with a time-varying
 * writer may not be removed.
 *
 * \param[in]  mesh_id  postprocessing mesh id
 */
/*----------------------------------------------------------------------------*/

void
cs_post_free_mesh(int  mesh_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check for the existence of a writer of the given id.
 *
 * \param[in]  writer_id  writer id to check
 *
 * \return  true if writer with this id exists, false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_post_writer_exists(int  writer_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check for the existence of a post-processing mesh of the given id.
 *
 * \param[in]  mesh_id  mesh id to check
 *
 * \return  true if mesh with this id exists, false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_post_mesh_exists(int  mesh_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the default writer format name
 *
 * \return  name of the default writer format
 */
/*----------------------------------------------------------------------------*/

const char *
cs_post_get_default_format(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the default writer format options
 *
 * \return  default writer format options string
 */
/*----------------------------------------------------------------------------*/

const char *
cs_post_get_default_format_options(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the next "reservable" (i.e. non-user) writer id available.
 *
 * \return  the smallest negative integer present, -1
 */
/*----------------------------------------------------------------------------*/

int
cs_post_get_free_writer_id(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the next "reservable" (i.e. non-user) mesh id available.
 *
 * \return  the smallest negative integer present, -1
 */
/*----------------------------------------------------------------------------*/

int
cs_post_get_free_mesh_id(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update "active" or "inactive" flag of writers based on the time step.
 *
 * Writers are activated if their output frequency is a divisor of the
 * current time step, or if their optional time step and value output lists
 * contain matches for the current time step.
 *
 * \param[in]  ts  time step status structure
 */
/*----------------------------------------------------------------------------*/

void
cs_post_activate_by_time_step(const cs_time_step_t  *ts);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Force the "active" or "inactive" flag for a specific writer or for all
 * writers for the current time step.
 *
 * \param[in]  writer_id  writer id, or 0 for all writers
 * \param[in]  activate   false to deactivate, true to activate
 */
/*----------------------------------------------------------------------------*/

void
cs_post_activate_writer(int   writer_id,
                        bool  activate);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Disable specific writer or all writers not currently active until
 *        \ref cs_post_enable_writer or \ref cs_post_activate_writer
 *        is called for those writers.
 *
 * For each call to this function for a given writer, the same number
 * of calls to \ref cs_post_enable_writer or a single call to
 * \ref cs_post_activate_writer is required to re-enable the writer.
 *
 * This is useful to disable output even of fixed meshes in preprocessing
 * stages.
 *
 * \param[in]  writer_id  writer id, or 0 for all writers
 */
/*----------------------------------------------------------------------------*/

void
cs_post_disable_writer(int   writer_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Enable a specific writer or all writers currently disabled by
 *        previous calls to \ref cs_post_disable_writer.
 *
 * For each previous call to \ref cs_post_disable_writer for a given writer,
 * a call to this function (or a single call to \ref cs_post_activate_writer)
 * is required to re-enable the writer.
 *
 * This is useful to disable output even of fixed meshes in preprocessing
 * stages.
 *
 * \param[in]  writer_id  writer id, or 0 for all writers
 */
/*----------------------------------------------------------------------------*/

void
cs_post_enable_writer(int   writer_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a pointer to the FVM writer associated to a writer_id.
 *
 * \param[in]  writer_id  associated writer id
 *
 * \return  a pointer to a fvm_writer_t structure
 */
/*----------------------------------------------------------------------------*/

fvm_writer_t *
cs_post_get_writer(int  writer_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return time dependency associated to a writer_id.
 *
 * \param[in]  writer_id  associated writer id
 *
 * \return  associated writer's time dependency
 */
/*----------------------------------------------------------------------------*/

fvm_writer_time_dep_t
cs_post_get_writer_time_dep(int  writer_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add an activation time step for a specific writer or for all writers.
 *
 * If a negative value is provided, a previously added activation time
 * step matching that absolute value will be removed, if present.
 *
 * \param[in]  writer_id  writer id, or 0 for all writers
 * \param[in]  nt         time step value to add (or remove)
 */
/*----------------------------------------------------------------------------*/

void
cs_post_add_writer_t_step(int  writer_id,
                          int  nt);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add an activation time value for a specific writer or for all writers.
 *
 * If a negative value is provided, a previously added activation time
 * step matching that absolute value will be removed, if present.
 *
 * \param[in]  writer_id  writer id, or 0 for all writers
 * \param[in]  t          time value to add (or remove)
 */
/*----------------------------------------------------------------------------*/

void
cs_post_add_writer_t_value(int     writer_id,
                           double  t);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Output post-processing meshes using associated writers.
 *
 * If the time step structure argument passed is NULL, a time-independent
 * output will be assumed.
 *
 * \param[in]  ts  time step status structure, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_post_write_meshes(const cs_time_step_t  *ts);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Output a variable defined at cells or faces of a post-processing mesh
 * using associated writers.
 *
 * \param[in]  mesh_id      id of associated mesh
 * \param[in]  writer_id    id of specified associated writer,
 *                          or \ref CS_POST_WRITER_ALL_ASSOCIATED for all
 * \param[in]  var_name     name of variable to output
 * \param[in]  var_dim      1 for scalar, 3 for vector, 6 for symmetric tensor,
 *                          9 for non-symmetric tensor
 * \param[in]  interlace    if a vector, true for interlaced values,
 *                          false otherwise
 * \param[in]  use_parent   true if values are defined on "parent" mesh,
 *                          false if values are defined on post-processing mesh
 * \param[in]  var_type     variable's data type
 * \param[in]  cel_vals     cell values
 * \param[in]  i_face_vals  interior face values
 * \param[in]  b_face_vals  boundary face values
 * \param[in]  ts           time step status structure, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_post_write_var(int                    mesh_id,
                  int                    writer_id,
                  const char            *var_name,
                  int                    var_dim,
                  bool                   interlace,
                  bool                   use_parent,
                  cs_post_type_t         var_type,
                  const void            *cel_vals,
                  const void            *i_face_vals,
                  const void            *b_face_vals,
                  const cs_time_step_t  *ts);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Output a variable defined at vertices of a post-processing mesh using
 *        associated writers.
 *
 * \param[in]  mesh_id     id of associated mesh
 * \param[in]  writer_id   id of specified associated writer,
 *                         or \ref CS_POST_WRITER_ALL_ASSOCIATED for all
 * \param[in]  var_name    name of variable to output
 * \param[in]  var_dim     1 for scalar, 3 for vector, 6 for symmetric tensor,
 *                         9 for non-symmetric tensor
 * \param[in]  interlace   if a vector, true for interlaced values,
 *                         false otherwise
 * \param[in]  use_parent  true if values are defined on "parent" mesh,
 *                         false if values are defined on post-processing mesh
 * \param[in]  var_type    variable's data type
 * \param[in]  vtx_vals    vertex values
 * \param[in]  ts          time step status structure, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_post_write_vertex_var(int                    mesh_id,
                         int                    writer_id,
                         const char            *var_name,
                         int                    var_dim,
                         bool                   interlace,
                         bool                   use_parent,
                         cs_post_type_t         var_type,
                         const void            *vtx_vals,
                         const cs_time_step_t  *ts);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Output an existing lagrangian particle attribute at particle
 *        positions or trajectory endpoints of a particle mesh using
 *        associated writers.
 *
 * \param[in]  mesh_id       id of associated mesh
 * \param[in]  writer_id     id of specified associated writer,
 *                           or \ref CS_POST_WRITER_ALL_ASSOCIATED for all
 * \param[in]  attr_id       associated particle attribute id
 * \param[in]  var_name      name of variable to output
 * \param[in]  component_id  if -1 : extract the whole attribute
 *                           if >0 : id of the component to extract
 * \param[in]  ts            time step status structure, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_post_write_particle_values(int                    mesh_id,
                              int                    writer_id,
                              int                    attr_id,
                              const char            *var_name,
                              int                    component_id,
                              const cs_time_step_t  *ts);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Output a variable defined at cells or faces of a post-processing mesh
 *        using associated writers.
 *
 * \param[in]  mesh_id              id of associated mesh
 * \param[in]  writer_id            id of specified associated writer,
 *                                  or \ref CS_POST_WRITER_ALL_ASSOCIATED for all
 * \param[in]  var_name             name of variable to output
 * \param[in]  var_dim              1 for scalar, 3 for vector, 6 for symmetric
 *                                  tensor, 9 for non-symmetric tensor
 * \param[in]  var_type             variable's data type
 * \param[in]  parent_location_id   asociated values mesh location, or 0
 *                                  if values are passed directly
 * \param[in]  interpolate_func     pointer to interpolation function,
 *                                  or NULL for default
 * \param[in]  interpolate_input    pointer to optional interpolation input
 *                                  data, or NULL for default
 * \param[in]  vals                 variable's values
 * \param[in]  ts                   time step status structure, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_post_write_probe_values(int                              mesh_id,
                           int                              writer_id,
                           const char                      *var_name,
                           int                              var_dim,
                           cs_post_type_t                   var_type,
                           int                              parent_location_id,
                           cs_interpolate_from_location_t  *interpolate_func,
                           void                            *interpolate_input,
                           const void                      *vals,
                           const cs_time_step_t            *ts);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update references to parent mesh of post-processing meshes in case of
 * computational mesh cell renumbering.
 *
 * This function may be called only once, after possible renumbering of cells,
 * to update existing post-processing meshes. Post-processing meshes defined
 * after renumbering will automatically be based upon the new numbering,
 * so this function will not need to be called again.
 *
 * \param[in]  init_cell_num  initial cell numbering (new -> old)
 */
/*----------------------------------------------------------------------------*/

void
cs_post_renum_cells(const cs_lnum_t  init_cell_num[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update references to parent mesh of post-processing meshes in case of
 * computational mesh interior and/or boundary faces renumbering.
 *
 * This function may be called only once, after possible renumbering of faces,
 * to update existing post-processing meshes. Post-processing meshes defined
 * after renumbering will automatically be based upon the new numbering,
 * so this function will not need to be called again.
 *
 * \param[in]  init_i_face_num  initial interior numbering (new -> old)
 * \param[in]  init_b_face_num  initial boundary numbering (new -> old)
 */
/*----------------------------------------------------------------------------*/

void
cs_post_renum_faces(const cs_lnum_t  init_i_face_num[],
                    const cs_lnum_t  init_b_face_num[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Configure the post-processing output so that mesh connectivity
 * may be automatically updated.
 *
 * This is done for meshes defined using selection criteria or functions.
 * The behavior of Lagrangian meshes is unchanged.
 *
 * To be effective, this function should be called before defining
 * postprocessing meshes.
 */
/*----------------------------------------------------------------------------*/

void
cs_post_set_changing_connectivity(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize post-processing writers
 */
/*----------------------------------------------------------------------------*/

void
cs_post_init_writers(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize main post-processing meshes
 *
 * The check_flag variable is a mask, used for additionnal post-processing:
 *
 *  - If (check_flag & 1), volume submeshes are output by groups if more
 *    than one group is present and the default writer uses the EnSight format.
 *
 *  - If (check_flag & 2), boundary submeshes are output by groups if more
 *    than one group is present and the default writer uses the EnSight format.
 *
 * It is recommended that post-processing meshes be defined before calling
 * this function, though specific "automatic" meshes (for example those
 * related to couplings) may be defined between this call and a time loop.
 *
 * \param[in]  check_mask  mask used for additional output
 */
/*----------------------------------------------------------------------------*/

void
cs_post_init_meshes(int check_mask);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check if post-processing is activated and then update post-processing
 *        of meshes if there is a need to update time-dependent meshes
 *
 * \param[in]  ts  time step status structure, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_post_time_step_begin(const cs_time_step_t  *ts);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Loop on post-processing meshes to output variables.
 *
 * This handles all default fields output, as well as all
 * registered output functions and outputs defined in
 * \ref cs_user_postprocess_values
 *
 * \param[in]  ts  time step status structure, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_post_time_step_output(const cs_time_step_t  *ts);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Flush writers and free time-varying and Lagragian mesh if needed
 *        of meshes if there is a time-dependent mesh
 */
/*----------------------------------------------------------------------------*/

void
cs_post_time_step_end(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Loop on post-processing meshes to output variables.
 *
 * This handles all default fields output, as well as all
 * registered output functions and outputs defined in
 * \ref cs_user_postprocess_values
 *
 * \param[in]  ts  time step status structure, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_post_write_vars(const cs_time_step_t  *ts);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy all structures associated with post-processing
 */
/*----------------------------------------------------------------------------*/

void
cs_post_finalize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Postprocess free (isolated) faces of the current global mesh
 */
/*----------------------------------------------------------------------------*/

void
cs_post_add_free_faces(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize post-processing writer with same format and associated
 * options as default writer, but no time dependency, intended to
 * troubleshoot errors.
 */
/*----------------------------------------------------------------------------*/

void
cs_post_init_error_writer(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize post-processing writer with same format and associated
 * options as default writer, but no time dependency, and associate
 * and output global volume mesh.
 *
 * This is intended to help troubleshoot errors using fields based
 * on cells.
 *
 * \return  id of error output mesh (< 0), or 0 if all writers are deactivated
 */
/*----------------------------------------------------------------------------*/

int
cs_post_init_error_writer_cells(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Register a processing of time-dependent variables to the call to
 * cs_post_write_vars().
 *
 * Note: if the input pointer is non-NULL, it must point to valid data
 * when the output function is called, so either:
 * - that value or structure should not be temporary (i.e. local);
 * - post-processing output must be ensured using cs_post_write_var()
 *   or similar before the data pointed to goes out of scope.
 *
 * \param[in]       function  function to register
 * \param[in, out]  input     pointer to optional (untyped) value or structure
 */
/*----------------------------------------------------------------------------*/

void
cs_post_add_time_dep_output(cs_post_time_dep_output_t  *function,
                            void                       *input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Register a processing of time-dependent variables than can be output
 * on different meshes to the call to cs_post_write_vars().
 *
 * Note: if the input pointer is non-NULL, it must point to valid data
 * when the output function is called, so either:
 * - that value or structure should not be temporary (i.e. local);
 * - post-processing output must be ensured using cs_post_write_var()
 *   or similar before the data pointed to goes out of scope.
 *
 * \param[in]       function  function to register
 * \param[in, out]  input     pointer to optional (untyped) value or structure
 */
/*----------------------------------------------------------------------------*/

void
cs_post_add_time_mesh_dep_output(cs_post_time_mesh_dep_output_t  *function,
                                 void                            *input);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_POST_H__ */
