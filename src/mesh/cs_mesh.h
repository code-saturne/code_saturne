#ifndef __CS_MESH_H__
#define __CS_MESH_H__

/*============================================================================
 * Main structure associated to a mesh
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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
#include "cs_halo.h"
#include "cs_interface.h"
#include "cs_numbering.h"
#include "cs_parall.h"
#include "cs_range_set.h"

#include "cs_mesh_builder.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*
 * Mesh modification type flags
 */

/*! Any type of mesh modification */
#define CS_MESH_MODIFIED (1 << 0)

/*! Mesh modification has changed mesh distribution balance */
#define CS_MESH_MODIFIED_BALANCE (1 << 1)

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! Mesh time dependency */
/*  -------------------- */

typedef enum {

  CS_MESH_FIXED,              /*!< Mesh definitions do not change with time */
  CS_MESH_TRANSIENT_COORDS,   /*!< Vertex coordinates may change with time */
  CS_MESH_TRANSIENT_CONNECT   /*!< Mesh connectivity may change with time */

} cs_mesh_time_dep_t;

/*! Mesh structure definition */
/*  ------------------------- */

typedef struct {

  /* General features */

  cs_lnum_t  dim;                /*!< space dimension */
  int        domain_num;         /*!< local domain number */
  int        n_domains;          /*!< number of domains */

  cs_mesh_time_dep_t  time_dep;  /*!< time dependency */

  /* Local dimensions */

  cs_lnum_t  n_cells;            /*!< number of cells */
  cs_lnum_t  n_i_faces;          /*!< number of interior faces */
  cs_lnum_t  n_b_faces;          /*!< number of boundary faces */
  cs_lnum_t  n_vertices;         /*!< number of vertices */

  cs_lnum_t  i_face_vtx_connect_size;  /*!< interior faces -> vertices
                                         connectivity size */
  cs_lnum_t  b_face_vtx_connect_size;  /*!< boundary faces -> vertices
                                         connectivity size */

  /* Local structures */

  cs_real_t    *vtx_coord;       /*!< vertex coordinates */

  cs_lnum_2_t  *i_face_cells;    /*!< interior faces -> cells connectivity */
  cs_lnum_t    *b_face_cells;    /*!< boundary faces -> cells connectivity */

  cs_lnum_t    *i_face_vtx_idx;  /*!< interior faces -> vertices index */
  cs_lnum_t    *i_face_vtx_lst;  /*!< interior faces -> vertices connectivity */

  cs_lnum_t    *b_face_vtx_idx;  /*!< boundary faces -> vertices index */
  cs_lnum_t    *b_face_vtx_lst;  /*!< boundary faces -> vertices connectivity */

  /* Global dimension */

  cs_gnum_t   n_g_cells;         /*!< global number of cells */
  cs_gnum_t   n_g_i_faces;       /*!< global number of interior faces */
  cs_gnum_t   n_g_b_faces;       /*!< global number of boundary faces */
  cs_gnum_t   n_g_vertices;      /*!< global number of vertices */

  cs_gnum_t   n_g_i_c_faces;     /*!< global number of interior faces
                                   for counts (with periodic faces
                                   counted only once) */

  /* Global numbering */

  cs_gnum_t  *global_cell_num;    /*!< global cell numbering */
  cs_gnum_t  *global_i_face_num;  /*!< global interior face numbering */
  cs_gnum_t  *global_b_face_num;  /*!< global boundary face numbering */
  cs_gnum_t  *global_vtx_num;     /*!< global vertex numbering */

  /* Periodictity features */

  int       n_init_perio;         /*!< number of initial periodicities */
  int       n_transforms;         /*!< number of transformations */

  int       have_rotation_perio;  /*!< periodicity rotation indicator */

  fvm_periodicity_t  *periodicity; /*!< parameters of each periodicity */

  /* Parallelism and/or periodic features */

  cs_halo_type_t  halo_type;       /*!< halo type */

  cs_lnum_t  n_cells_with_ghosts;  /*!< total number of cells on the local rank
                                        (n_cells + n_ghost_cells) */
  cs_lnum_t  n_ghost_cells;        /*!< number of "ghost" cells */

  cs_interface_set_t  *vtx_interfaces;  /*!< vertices interface set */
  cs_halo_t           *halo;            /*!< ghost cells structure */
  cs_range_set_t      *vtx_range_set;   /*!< handle local/distant ranges for
                                          vertices in parallel */

  cs_numbering_t  *cell_numbering;      /*!< cell numbering info */
  cs_numbering_t  *vtx_numbering;       /*!< vertex numbering info */
  cs_numbering_t  *i_face_numbering;    /*!< interior face numbering info */
  cs_numbering_t  *b_face_numbering;    /*!< boundary face numbering info */

  /* Re-computable connectivity features */

  cs_lnum_t   n_b_cells;             /*!< number of boundary cells */
  cs_lnum_t  *b_cells;               /*!< boundary cell list */

  /* Extended neighborhood features */

  cs_lnum_t  *cell_cells_idx;  /*!< "cell -> cells" connectivity index for
                                 extended halo. Only defined if extended
                                 neighborhood is built. */
  cs_lnum_t  *cell_cells_lst;  /*!< "cell -> cells" connectivity list for
                                 extended halo. Only defined if extended
                                 neighborhood is built. */

  cs_lnum_t  *gcell_vtx_idx;   /*!< ghost cells -> vertices index */
  cs_lnum_t  *gcell_vtx_lst;   /*!< ghost cells -> vertices list */

  /* Group and family features */

  int         n_groups;            /*!< number of groups */
  int        *group_idx;           /*!< starting index in group */
  char       *group;               /*!< list of group names */

  int         n_families;          /*!< number of families */
  int         n_max_family_items;  /*!< max. number of items for one family */
  int        *family_item;         /*!< family items */
  int        *cell_family;         /*!< cell family */
  int        *i_face_family;       /*!< interior face family */
  int        *b_face_family;       /*!< boundary face family */

  fvm_group_class_set_t *class_defs;  /*!< definition of group classes for
                                        selection and postprocessing (built
                                        from element families and their
                                        descriptions) */
  fvm_selector_t  *select_cells;      /*!< cells selection object */
  fvm_selector_t  *select_i_faces;    /*!< interior faces selection object */
  fvm_selector_t  *select_b_faces;    /*!< boundary faces selection object */

  /* Refinement features */

  bool        have_r_gen;          /*!< has mesh refinement information */
  char       *i_face_r_gen;        /*!< interior face refinement generation */
  char       *vtx_r_gen;           /*!< vertex refinement generation */

  /* Status flags */

  cs_gnum_t n_g_free_faces;        /*!< global number of boundary faces
                                     which are in fact isolated */

  cs_gnum_t n_g_b_faces_all;       /*!< global number of boundary faces
                                     including those ignored in FV schemes */
  cs_lnum_t n_b_faces_all;         /*!< number of boundary faces including
                                     faces ignored in FV schemes */

  int verbosity;                   /*!< current verbosity level */
  int modified;                    /*!< modification status */
  int save_if_modified;            /*!< flag for mesh saving behavior:
                                     0: never save
                                     1: saved when modified (default)
                                     2: always save */

} cs_mesh_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

extern cs_mesh_t *cs_glob_mesh; /* Pointer to main mesh structure */

/*============================================================================
 *  Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Update a scalar array in case of parallelism and/or periodicity.
 *
 * Fortran interface:
 *
 * subroutine synsca(var)
 * *****************
 *
 * var   : <-> : scalar array
 *----------------------------------------------------------------------------*/

void CS_PROCF(synsca, SYNSCA)
(
 cs_real_t  var[]
);

/*----------------------------------------------------------------------------
 * Update a scalar array in case of parallelism and/or periodicity,
 * using an extended halo.
 *
 * Fortran interface:
 *
 * subroutine synsce(var)
 * *****************
 *
 * var   : <-> : scalar array
 *----------------------------------------------------------------------------*/

void CS_PROCF(synsce, SYNSCE)
(
 cs_real_t  var[]
);

/*----------------------------------------------------------------------------
 * Update a vector array in case of parallelism and/or periodicity.
 *
 * Fortran interface:
 *
 * subroutine synvin(var)
 * *****************
 *
 * var   : <-> : interleaved vector (of dimension 3)
 *----------------------------------------------------------------------------*/

void CS_PROCF(synvin, SYNVIN)
(
 cs_real_t  var[]
);

/*----------------------------------------------------------------------------
 * Update a vector array in case of parallelism and/or periodicity,
 * using an extended halo.
 *
 * Fortran interface:
 *
 * subroutine synvin(var)
 * *****************
 *
 * var   : <-> : interleaved vector (of dimension 3)
 *----------------------------------------------------------------------------*/

void CS_PROCF(synvie, SYNVIE)
(
 cs_real_t  var[]
);

/*----------------------------------------------------------------------------
 * Update a tensor array in case of parallelism and/or periodicity.
 *
 * Fortran interface:
 *
 * subroutine syntin(var)
 * *****************
 *
 * var   : <-> : interleaved tensor (of dimension 3x3)
 *----------------------------------------------------------------------------*/

void CS_PROCF(syntin, SYNTIN)
(
 cs_real_t  var[]
);

/*----------------------------------------------------------------------------
 * Update a symmetric tensor array in case of parallelism and/or periodicity.
 *
 * Fortran interface:
 *
 * subroutine syntis(var)
 * *****************
 *
 * var   : <-> : interleaved symmetric tensor (of dimension 6)
 *----------------------------------------------------------------------------*/

void CS_PROCF(syntis, SYNTIS)
(
 cs_real_t  var[]
);

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create an empty mesh structure
 *
 * returns:
 *   pointer to created mesh structure
 *----------------------------------------------------------------------------*/

cs_mesh_t *
cs_mesh_create(void);

/*----------------------------------------------------------------------------
 * Destroy a mesh structure
 *
 * mesh       <->  pointer to a mesh structure
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

cs_mesh_t *
cs_mesh_destroy(cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------
 * Reinitialize mesh structure.
 *
 * returns:
 *   pointer to created mesh structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_reinit(cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------
 * Update (compactify) an array of global numbers.
 *
 * parameters:
 *   n_elts   <-> number of local elements
 *   elt_gnum <-> global element numbers
 *
 * return:
 *   associated global number of elements
 *----------------------------------------------------------------------------*/

cs_gnum_t
cs_mesh_compact_gnum(cs_lnum_t   n_elts,
                     cs_gnum_t  *elt_gnum);

/*----------------------------------------------------------------------------
 * Remove arrays and structures that mey be rebuilt.
 *
 * mesh       <-> pointer to a mesh structure
 * free_halos <-- if true, free halos and parallel/periodic interface
 *                structures
 *----------------------------------------------------------------------------*/

void
cs_mesh_free_rebuildable(cs_mesh_t  *mesh,
                         bool        free_halos);

/*----------------------------------------------------------------------------
 * Discard free (isolated) faces from a mesh.
 *
 * This should always be done before using the mesh for computation.
 *
 * parameters:
 *   mesh  <->  pointer to mesh structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_discard_free_faces(cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------
 * Discard free (isolated) vertices from a mesh.
 *
 * This is recommended before using the mesh for computation.
 *
 * parameters:
 *   mesh    <-> pointer to mesh structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_discard_free_vertices(cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Discard mesh refinement info.
 *
 * This information is used only for mesh coarsening or post-processing output
 * of the refinement level, so can be discarded in other cases.
 *
 * \param[in, out]  mesh  pointer to mesh structure
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_discard_refinement_info(cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------
 * Generate or update list of mesh boundary cells.
 *
 * parameters:
 *   mesh  <->  pointer to a cs_mesh_t structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_update_b_cells(cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------
 * Compute or update mesh structure members that depend on other members,
 * but whose results may be reused, such as global number of elements
 * (cells, vertices, interior and boundary faces) and sync cell family.
 *
 * parameters:
 *   mesh   <->  pointer to a cs_mesh_t structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_update_auxiliary(cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------
 * Creation and initialization of mesh face and vertex interfaces.
 *
 * parameters:
 *   mesh  <->  pointer to mesh structure
 *   mb    <->  pointer to mesh builder (in case of periodicity)
 *----------------------------------------------------------------------------*/

void
cs_mesh_init_interfaces(cs_mesh_t          *mesh,
                        cs_mesh_builder_t  *mb);

/*----------------------------------------------------------------------------*/
/*!
 * Creation and initialization of halo structures.
 *
 * Treatment of parallel and/or periodic halos for standard and extended
 * ghost cells according to halo type requested by global options.
 *
 * \param[in, out]  mesh                   pointer to mesh structure
 * \param[in, out]  mb                     pointer to mesh builder
 *                                         (for periodicity)
 * \param[in]       halo_type              type of halo (standard or extended)
 * \param[in]       verbosity              verbosity
 * \param[in]       rebuild_vtx_interface  also rebuild vertex interfaces ?
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_init_halo(cs_mesh_t          *mesh,
                  cs_mesh_builder_t  *mb,
                  cs_halo_type_t      halo_type,
                  int                 verbosity,
                  bool                rebuild_vtx_interace);

/*----------------------------------------------------------------------------
 * Get the global number of ghost cells.
 *
 * parameters:
 *  mesh <--  pointer to a mesh structure
 *
 * returns:
 *  Global number of ghost cells
 *---------------------------------------------------------------------------*/

cs_gnum_t
cs_mesh_n_g_ghost_cells(cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------
 * Update a scalar array in case of parallelism and/or periodicity.
 *
 * Note: this function is only present so that a C equivalent to the
 *       Fortran wrappers is available. In C code, directly using
 *       cs_halo_sync_var() is preferred.
 *
 * parameters:
 *   var  <->  scalar array
 *----------------------------------------------------------------------------*/

void
cs_mesh_sync_var_scal(cs_real_t  *var);

/*----------------------------------------------------------------------------
 * Update a scalar array in case of parallelism and/or periodicity,
 * using an extended halo.
 *
 * Note: this function is only present so that a C equivalent to the
 *       Fortran wrappers is available. In C code, directly using the
 *       cs_halo_sync_var() is preferred.
 *
 * parameters:
 *   var  <->  scalar array
 *----------------------------------------------------------------------------*/

void
cs_mesh_sync_var_scal_ext(cs_real_t  *var);

/*----------------------------------------------------------------------------
 * Update a vector array in case of parallelism and/or periodicity.
 *
 * parameters:
 *   var  <->  interleaved vector (of dimension 3)
 *----------------------------------------------------------------------------*/

void
cs_mesh_sync_var_vect(cs_real_t  *var);

/*----------------------------------------------------------------------------
 * Update a vector array in case of parallelism and/or periodicity,
 * using an extended halo.
 *
 * parameters:
 *   var  <->  interleaved vector (of dimension 3)
 *----------------------------------------------------------------------------*/

void
cs_mesh_sync_var_vect_ext(cs_real_t  *var);

/*----------------------------------------------------------------------------
 * Update a tensor array in case of parallelism and/or periodicity.
 *
 * parameters:
 *   var  <->  interleaved tensor (of dimension 3x3)
 *----------------------------------------------------------------------------*/

void
cs_mesh_sync_var_tens(cs_real_t  *var);

/*----------------------------------------------------------------------------
 * Update a symmetric tensor array in case of parallelism and/or periodicity.
 *
 * parameters:
 *   var  <->  symmetric interleaved tensor (of dimension 6)
 *----------------------------------------------------------------------------*/

void
cs_mesh_sync_var_sym_tens(cs_real_6_t  *var);

/*----------------------------------------------------------------------------
 * Order family numbers and remove duplicates
 *
 * parameters
 *   mesh <-> pointer to mesh structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_clean_families(cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------
 * Create group classes based on a mesh's family definitions.
 *
 * parameters:
 *   mesh <-- pointer to mesh structure
 *
 * returns:
 *   pointer to group classes structure based on mesh's family definitions
 *----------------------------------------------------------------------------*/

fvm_group_class_set_t *
cs_mesh_create_group_classes(cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------
 * Define group classes for a mesh based on its family definitions.
 *
 * parameters:
 *   mesh <-> pointer to mesh structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_init_group_classes(cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------
 * Assign selectors to global mesh.
 *
 * Should be called once the mesh is fully built.
 *----------------------------------------------------------------------------*/

void
cs_mesh_init_selectors(void);

/*----------------------------------------------------------------------------
 * Update selector and associated structures.
 *
 * parameters:
 *   mesh  <-> pointer to a mesh structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_update_selectors(cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------
 * Get global lists of periodic face couples.
 *
 * In parallel, each face couple may appear on only one rank.
 *
 * The caller is responsible for freeing the arrays allocated and returned
 * by this function once they are no onger needed.
 *
 * parameters:
 *   mesh                 <-- pointer to mesh structure
 *   n_perio_face_couples --> global number of periodic couples per
 *                            periodicity (size: mesh->n_init_perio)
 *   perio_face_couples   --> arrays of global periodic couple face numbers,
 *                            for each periodicity
 *----------------------------------------------------------------------------*/

void
cs_mesh_get_perio_faces(const cs_mesh_t    *mesh,
                        cs_lnum_t         **n_perio_face_couples,
                        cs_gnum_t        ***perio_face_couples);

/*----------------------------------------------------------------------------
 * Build global cell numbering array extended to ghost cell values.
 *
 * If the blank_perio flag is nonzero, periodic ghost cell numbers
 * are set to zero instead of the value of the matching cell.
 *
 * The caller is responsible for freeing the returned array when it
 * is no longer useful.
 *
 * parameters:
 *   mesh        <-- pointer to mesh structure
 *   blank_perio <-- flag to zeroe periodic cell values
 *----------------------------------------------------------------------------*/

cs_gnum_t *
cs_mesh_get_cell_gnum(const cs_mesh_t  *mesh,
                      int               blank_perio);

/*----------------------------------------------------------------------------
 * Mark interior faces with the number of their associated periodic
 * transform id.
 *
 * parameters:
 *   mesh      <-- pointer to mesh structure
 *   perio_num --> periodicity number associated with each face, signed for
 *                 direct/reverse transform, 0 for non-periodic faces
 *                 (size: mesh->n_i_faces)
 *----------------------------------------------------------------------------*/

void
cs_mesh_get_face_perio_num(const cs_mesh_t  *mesh,
                           int               perio_num[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Mark cells adjacent to boundary, through faces or vertices.
 *
 * Note that cells adjacent through a boundary face can be accessed
 * directly through the mesh->b_cells member, but the set of cells flagged
 * by this function is more complete.
 *
 * \param[in]   mesh         pointer to a mesh structure
 * \param[out]  cell_b_flag  1 for cells adjacent to boundary, 0 for others
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_tag_boundary_cells(cs_mesh_t  *mesh,
                           int         cell_b_flag[]);

/*----------------------------------------------------------------------------
 * Print information on a mesh structure.
 *
 * parameters:
 *   mesh  <--  pointer to mesh structure.
 *   name  <--  associated name.
 *----------------------------------------------------------------------------*/

void
cs_mesh_print_info(const cs_mesh_t  *mesh,
                   const char       *name);

/*----------------------------------------------------------------------------
 * Compute global face connectivity size.
 *
 * Faces on simple parallel boundaries are counted only once, but periodic
 * faces are counted twice.
 *
 * parameters:
 *   mesh                   <-- pointer to a cs_mesh_t structure
 *   g_i_face_vertices_size --> global interior face connectivity size, or NULL
 *   g_b_face_vertices_size --> global boundary face connectivity size, or NULL
 *----------------------------------------------------------------------------*/

void
cs_mesh_g_face_vertices_sizes(const cs_mesh_t  *mesh,
                              cs_gnum_t        *g_i_face_vertices_size,
                              cs_gnum_t        *g_b_face_vertices_size);

/*----------------------------------------------------------------------------
 * Print statistics about mesh selectors usage to log.
 *
 * parameters:
 *   mesh <-- pointer to a mesh structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_selector_stats(cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------
 * Dump of a mesh structure.
 *
 * parameters:
 *   mesh  <->  pointer to mesh structure.
 *----------------------------------------------------------------------------*/

void
cs_mesh_dump(const cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Determine number of blocks and associated groups to be used
 *        for loops on interior faces.
 *
 * Blocks from a same group may be processed in parallel, while blocks from
 * separate groups must not run simultaneously to avoid race conditions.
 *
 * \param[in]   m           pointer to mesh
 * \param[in]   e2n         associated indexed sum algorithm type.
 * \param[in]   block_size  size of thread blocks (chunks) if > 0,
 *                          ignored (recomputed) if 0.
 * \param[out]  n_groups    number of associated block groups
 * \param[out]  n_blocks    number of associated blocks
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_i_faces_thread_block_count(const cs_mesh_t     *m,
                                   const cs_e2n_sum_t   e2n,
                                   int                  block_size,
                                   int                 *n_groups,
                                   int                 *n_blocks);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute array index bounds for a block of interior faces
 *        associated to a thread or task.
 *
 * When the CS_E2N_SUM_SCATTER indexed sum algorithmm is used and mesh
 * interior faces are renumbered for threads, the bounds provided are
 * those based on the matching group and thread.
 *
 * In other cases, if block_size < 1 (i.e. not specified), the start and
 * past-the-end indexes are defined so as to evenly distribute values
 * to threads, in a manner similar to OpenMP <tt>schedule(static)</tt>.
 * With a block size larger than zero, indexes will be simply based on
 * the block's start and past-the end index.
 *
 * \param[in]       m            pointer to mesh
 * \param[in]       e2n          associated indexed sum algorithm type.
 * \param[in]       group_id     group id
 * \param[in]       block_id     block id (relative to group)
 * \param[in]       block_count  number of blocks
 * \param[in]       block_size   size of blocks (chunks) if > 0,
 *                               ignored (recomputed) if 0
 * \param[in, out]  s_id         start index for the current block
 * \param[in, out]  e_id         past-the-end index for the current block
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_i_faces_thread_block_range(const cs_mesh_t     *m,
                                   const cs_e2n_sum_t   e2n,
                                   int                  group_id,
                                   int                  block_id,
                                   int                  block_count,
                                   int                  block_size,
                                   cs_lnum_t           *s_id,
                                   cs_lnum_t           *e_id);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MESH_H__ */
