#ifndef __CS_MESH_H__
#define __CS_MESH_H__

/*============================================================================
 * Main structure associated to a mesh
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
#include "cs_range_set.h"

#include "cs_mesh_builder.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Mesh structure definition */
/* ------------------------- */

typedef struct {

  /* General features */

  cs_lnum_t  dim;                  /* Space dimension */
  cs_lnum_t  domain_num;           /* Local domain number */
  cs_lnum_t  n_domains;            /* Number of domains */

  /* Local dimensions */

  cs_lnum_t  n_cells;              /* Number of cells */
  cs_lnum_t  n_i_faces;            /* Number of interior faces */
  cs_lnum_t  n_b_faces;            /* Number of boundary faces */
  cs_lnum_t  n_vertices;           /* Number of vertices */

  cs_lnum_t  i_face_vtx_connect_size;  /* Size of the connectivity
                                          interior faces -> vertices */
  cs_lnum_t  b_face_vtx_connect_size;  /* Size of the connectivity
                                          boundary faces -> vertices */

  /* Local structures */

  cs_real_t    *vtx_coord;         /* Vertex coordinates */

  cs_lnum_2_t  *i_face_cells;      /* Interior faces -> cells connectivity */
  cs_lnum_t    *b_face_cells;      /* Boundary faces -> cells connectivity */

  cs_lnum_t    *i_face_vtx_idx;    /* Interior faces -> vertices index */
  cs_lnum_t    *i_face_vtx_lst;    /* Interior faces -> vertices connectivity */

  cs_lnum_t    *b_face_vtx_idx;    /* Boundary faces -> vertices index */
  cs_lnum_t    *b_face_vtx_lst;    /* Boundary faces -> vertices connectivity */

  /* Global dimension */

  cs_gnum_t   n_g_cells;           /* Global number of cells */
  cs_gnum_t   n_g_i_faces;         /* Global number of interior faces */
  cs_gnum_t   n_g_b_faces;         /* Global number of boundary faces */
  cs_gnum_t   n_g_vertices;        /* Global number of vertices */

  /* Global numbering */

  cs_gnum_t  *global_cell_num;     /* Global cell numbering */
  cs_gnum_t  *global_i_face_num;   /* Global interior face numbering */
  cs_gnum_t  *global_b_face_num;   /* Global boundary face numbering */
  cs_gnum_t  *global_vtx_num;      /* Global vertex numbering */

  /* Periodictity features */

  int       n_init_perio;          /* Number of initial periodicities */
  int       n_transforms;          /* Number of transformations */

  int       have_rotation_perio;   /* Periodicity rotation indicator */

  fvm_periodicity_t  *periodicity; /* parameters of each periodicity */

  /* Parallelism and/or periodic features */

  cs_halo_type_t  halo_type;         /* Halo type */

  cs_lnum_t  n_cells_with_ghosts;    /* Total number of cells on the local rank
                                        (n_cells + n_ghost_cells) */
  cs_lnum_t  n_ghost_cells;          /* Number of "ghost" cells */

  cs_interface_set_t  *vtx_interfaces;   /* Vertices interface set */
  cs_halo_t           *halo;             /* Ghost cells structure */
  cs_range_set_t      *vtx_range_set;    /* Handle local/distant ranges for
                                            vertices in parallel */

  cs_numbering_t  *cell_numbering;   /* Cell numbering info */
  cs_numbering_t  *i_face_numbering; /* Interior face numbering info */
  cs_numbering_t  *b_face_numbering; /* Boundary face numbering info */

  /* Re-computable connectivity features */

  cs_lnum_t   n_b_cells;             /* Number of boundary cells */
  cs_lnum_t  *b_cells;               /* Boundary cell list */

  /* Extended neighborhood features */

  cs_lnum_t  *cell_cells_idx;  /* "cell -> cells" connectivity index for
                                  extended halo. Only defined if extended
                                  neighborhood is built. */
  cs_lnum_t  *cell_cells_lst;  /* "cell -> cells" connectivity list for
                                  extended halo. Only defined if extended
                                  neighborhood is built. */

  cs_lnum_t  *gcell_vtx_idx;   /* ghost cells -> vertices index */
  cs_lnum_t  *gcell_vtx_lst;   /* ghost cells -> vertices list */

  /* Group and family features */

  int         n_groups;            /* Number of groups */
  int        *group_idx;           /* Starting index in group */
  char       *group;               /* List of group names */

  int         n_families;          /* Number of families */
  int         n_max_family_items;  /* Max. number of items for one family */
  int        *family_item;         /* Family items */
  int        *cell_family;         /* Cell family */
  int        *i_face_family;       /* Interior face family */
  int        *b_face_family;       /* Boundary face family */

  fvm_group_class_set_t *class_defs;  /* Definition of group classes for
                                         selection and postprocessing (built
                                         from element families and their
                                         descriptions) */
  fvm_selector_t  *select_cells;      /* Cells selection object */
  fvm_selector_t  *select_i_faces;    /* Interior faces selection object */
  fvm_selector_t  *select_b_faces;    /* Boundary faces selection object */

  /* Status flags */

  cs_gnum_t n_g_free_faces;          /* Global number of boundary faces
                                        which are in fact isolated */
  int verbosity;                     /* Current verbosity level */
  int modified;                      /* Modification status */

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
 * Update a scalar array in case of parallelism and/or periodicity,
 * ignoring periodicity of rotation
 *
 * Fortran interface:
 *
 * subroutine syncmp(var)
 * *****************
 *
 * var   : <-> : scalar array
 *----------------------------------------------------------------------------*/

void CS_PROCF(syncmp, SYNCMP)
(
 cs_real_t  var[]
);

/*----------------------------------------------------------------------------
 * Update a vector array in case of parallelism and/or periodicity.
 *
 * Fortran interface:
 *
 * subroutine synvec(var)
 * *****************
 *
 * var1   : <-> : vector component 1 array
 * var2   : <-> : vector component 2 array
 * var3   : <-> : vector component 3 array
 *----------------------------------------------------------------------------*/

void CS_PROCF(synvec, SYNVEC)
(
 cs_real_t  var1[],
 cs_real_t  var2[],
 cs_real_t  var3[]
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
 * Update a diagonal tensor array in case of parallelism and/or periodicity.
 *
 * Fortran interface:
 *
 * subroutine syndia(var)
 * *****************
 *
 * var11   : <-> : diagonal tensor component 11 array
 * var22   : <-> : diagonal tensor component 22 array
 * var33   : <-> : diagonal tensor component 33 array
 *----------------------------------------------------------------------------*/

void CS_PROCF(syndia, SYNDIA)
(
 cs_real_t  var11[],
 cs_real_t  var22[],
 cs_real_t  var33[]
);

/*----------------------------------------------------------------------------
 * Update a tensor array in case of parallelism and/or periodicity.
 *
 * Fortran interface:
 *
 * subroutine synten(var)
 * *****************
 *
 * var11   : <-> : tensor component 11 array
 * var12   : <-> : tensor component 12 array
 * var13   : <-> : tensor component 13 array
 * var21   : <-> : tensor component 21 array
 * var22   : <-> : tensor component 22 array
 * var23   : <-> : tensor component 23 array
 * var31   : <-> : tensor component 31 array
 * var32   : <-> : tensor component 32 array
 * var33   : <-> : tensor component 33 array
 *----------------------------------------------------------------------------*/

void CS_PROCF(synten, SYNTEN)
(
 cs_real_t  var11[],
 cs_real_t  var12[],
 cs_real_t  var13[],
 cs_real_t  var21[],
 cs_real_t  var22[],
 cs_real_t  var23[],
 cs_real_t  var31[],
 cs_real_t  var32[],
 cs_real_t  var33[]
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

/*----------------------------------------------------------------------------
 * Creation and initialization of halo structures.
 *
 * Treatment of parallel and/or periodic halos for standard and extended
 * ghost cells according to halo type requested by global options.
 *
 * parameters:
 *   mesh       <->  pointer to mesh structure
 *   mb         <->  pointer to mesh builder (in case of periodicity)
 *   halo_type  <->  type of halo (standard or extended)
 *----------------------------------------------------------------------------*/

void
cs_mesh_init_halo(cs_mesh_t          *mesh,
                  cs_mesh_builder_t  *mb,
                  cs_halo_type_t      halo_type);

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
 * Update a component of a vector for parallelism and/or periodicity,
 * ignoring periodicity of rotation.
 *
 * Note: this function is only present so that a C equivalent to the
 *       Fortran wrappers is available. In C code, directly using the
 *       cs_halo_sync_var() is preferred.
 *
 * parameters:
 *   var  <->  scalar array
 *----------------------------------------------------------------------------*/

void
cs_mesh_sync_var_component(cs_real_t  *var);

/*----------------------------------------------------------------------------
 * Update a vector array in case of parallelism and/or periodicity.
 *
 * parameters:
 *   var1  <->  vector component 1 array
 *   var2  <->  vector component 2 array
 *   var3  <->  vector component 3 array
 *----------------------------------------------------------------------------*/

void
cs_mesh_sync_var_vect_ni(cs_real_t  *var1,
                         cs_real_t  *var2,
                         cs_real_t  *var3);

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
 * Update a components of a vector for parallelism and/or periodicity,
 * ignoring periodicity of rotation.
 *
 *   var                  <-> gradient components (interleaved)
 *----------------------------------------------------------------------------*/

void
cs_mesh_sync_var_vect_no_rotation(cs_real_t  *var);

/*----------------------------------------------------------------------------
 * Update a diagonal tensor array in case of parallelism and/or periodicity.
 *
 * parameters:
 *   var11  <->  diagonal tensor component 11 array
 *   var22  <->  diagonal tensor component 22 array
 *   var33  <->  diagonal tensor component 33 array
 *----------------------------------------------------------------------------*/

void
cs_mesh_sync_var_diag_ni(cs_real_t  *var11,
                         cs_real_t  *var22,
                         cs_real_t  *var33);

/*----------------------------------------------------------------------------
 * Update a tensor array in case of parallelism and/or periodicity.
 *
 * parameters:
 *   var11  <->  tensor component 11 array
 *   var12  <->  tensor component 12 array
 *   var13  <->  tensor component 13 array
 *   var21  <->  tensor component 21 array
 *   var22  <->  tensor component 22 array
 *   var23  <->  tensor component 23 array
 *   var31  <->  tensor component 31 array
 *   var32  <->  tensor component 32 array
 *   var33  <->  tensor component 33 array
 *----------------------------------------------------------------------------*/

void
cs_mesh_sync_var_tens_ni(cs_real_t  *var11,
                         cs_real_t  *var12,
                         cs_real_t  *var13,
                         cs_real_t  *var21,
                         cs_real_t  *var22,
                         cs_real_t  *var23,
                         cs_real_t  *var31,
                         cs_real_t  *var32,
                         cs_real_t  *var33);

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
cs_mesh_sync_var_sym_tens(cs_real_t  *var);

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

END_C_DECLS

#endif /* __CS_MESH_H__ */
