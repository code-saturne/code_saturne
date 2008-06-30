/*============================================================================
 *
 *                    Code_Saturne version 1.3
 *                    ------------------------
 *
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2008 EDF S.A., France
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

#ifndef __CS_MESH_H__
#define __CS_MESH_H__

/*============================================================================
 * Main structure associated to a mesh
 *============================================================================*/

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <fvm_defs.h>
#include <fvm_group.h>
#include <fvm_selector.h>
#include <fvm_periodicity.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/* Halo type */

typedef enum {

  CS_MESH_HALO_STANDARD,
  CS_MESH_HALO_EXTENDED,
  CS_MESH_HALO_N_TYPES

} cs_mesh_halo_type_t;

/*============================================================================
 * Type definition
 *============================================================================*/

/* Structure for halo management */
/* ----------------------------- */

typedef struct {

  cs_int_t  n_c_domains;     /* Number of communicating domains. */
  cs_int_t  *c_domain_rank;  /* List of communicating ranks */

  /* in_halo features : send to distant ranks */

  cs_int_t  n_elts_in[2];    /* Numer of ghost elements in in_halo
                                n_elts[0] = standard elements
                                n_elts[1] = extended + standard elements */

  cs_int_t  *list_in;        /* List of local numbers of elements in in_halo */

  cs_int_t  *index_in;       /* Index on in_elements.
                                Size = 2*n_c_domains. For each rank, we
                                have an index for standard halo and one
                                for extended halo. */

  cs_int_t  *perio_lst_in;  /* For each transformation and for each type of halo
                               on each communicating rank, we store 2 data:
                                 - start index,
                                 - number of elements. */

  /* out_halo features : receive from distant ranks */

  cs_int_t  n_elts_out[2];  /* Numer of ghost elements in out_halo
                                 n_elts[0] = standard elements
                                 n_elts[1] = extended + standard elements */

  cs_int_t  *list_out;      /* List of local numbers of elements in out_halo */

  cs_int_t  *index_out;     /* Index on in_elements.
                               Size = 2*n_c_domains. For each rank, we
                               have an index for standard halo and one
                               for extended halo. */

  cs_int_t  *perio_lst_out; /* For each transformation and for each type of halo
                               on each communicating rank, we store 2 data:
                                 - start index,
                                 - number of elements. */

  cs_real_t  *tmp_buffer;   /* Buffer used to de-interlace variable
                               in case of strided variable to sync. */

#if defined(_CS_HAVE_MPI)
  MPI_Request   *mpi_request;   /* MPI Request array */
  MPI_Status    *mpi_status;    /* MPI Status array */

  cs_real_t  *comm_buffer;      /* Buffer for the communication purpose.
                                   Buffer size is equal to the maximum
                                   number of ghost cells between in_halo and
                                   out_halo. */
#endif

  /* Organisation of perio_lst:

         -------------------------------------------------
    T1:  |   |   |   |   |   |   |   |   |   |   |   |   |
         -------------------------------------------------
          idx  n  idx  n  idx  n  idx  n  idx  n  idx  n
          ______  ______  ______  ______  ______  ______
           std     ext     std     ext     std     ext
           ___________     ___________     ___________
             rank 0          rank 1          rank 2

         -------------------------------------------------
    T2:  |   |   |   |   |   |   |   |   |   |   |   |   |
         -------------------------------------------------
          idx  n  idx  n  idx  n  idx  n  idx  n  idx  n
          ______  ______  ______  ______  ______  ______
           std     ext     std     ext     std     ext
           ___________     ___________     ___________
             rank 0          rank 1          rank 2

         -------------------------------------------------
    T3:  |   |   |   |   |   |   |   |   |   |   |   |   |
         -------------------------------------------------
          idx  n  idx  n  idx  n  idx  n  idx  n  idx  n
          ______  ______  ______  ______  ______  ______
           std     ext     std     ext     std     ext
           ___________     ___________     ___________
             rank 0          rank 1          rank 2

  etc...

  */

} cs_mesh_halo_t;


/* Mesh structure definition */
/* ------------------------- */

typedef struct {

  /* General features */

  cs_int_t   dim;                  /* Space dimension */
  cs_int_t   domain_num;           /* Local domain number */
  cs_int_t   n_domains;            /* Number of domains */


  /* Local dimensions */

  cs_int_t   n_cells;              /* Number of cells */
  cs_int_t   n_i_faces;            /* Number of internal faces */
  cs_int_t   n_b_faces;            /* Number of border faces */
  cs_int_t   n_vertices;           /* Number of vertices */

  cs_int_t   i_face_vtx_connect_size;  /* Size of the connectivity
                                          internal faces -> vertices */
  cs_int_t   b_face_vtx_connect_size;  /* Size of the connectivity
                                          border faces -> vertices */

  /* Local structures */

  cs_real_t  *vtx_coord;          /* Coordinates of the vertices */

  cs_int_t   *i_face_cells;       /* Internal faces -> cells connectivity */
  cs_int_t   *b_face_cells;       /* Border faces -> cells connectivity */

  cs_int_t   *i_face_vtx_idx;     /* Internal faces -> vertices connectivity
                                     index */
  cs_int_t   *i_face_vtx_lst;     /* Internal faces -> vertices connectivity */

  cs_int_t   *b_face_vtx_idx;     /* Border faces -> vertices connectivity
                                     index */
  cs_int_t   *b_face_vtx_lst;     /* Border faces -> vertices connectivity */


  /* Global dimension */

  fvm_gnum_t   n_g_cells;           /* Global number of cells */
  fvm_gnum_t   n_g_i_faces;         /* Global number of internal faces */
  fvm_gnum_t   n_g_b_faces;         /* Global number of border faces */
  fvm_gnum_t   n_g_vertices;        /* Global number of vertices */

  /* Global numbering */

  fvm_gnum_t  *global_cell_num;    /* Global cell numbering */
  fvm_gnum_t  *global_i_face_num;  /* Global internal face numbering
                                      (remain in the initial local numbering) */
  fvm_gnum_t  *global_b_face_num;  /* Global border face numbering
                                      (remain in the initial local numbering) */
  fvm_gnum_t  *global_vtx_num;     /* Global vertex numbering */


  /* Initial face numbering (for exchange and restart) */

  cs_int_t   *init_i_face_num;    /* Initial internal face numbering
                                    (if renumbering) */
  cs_int_t   *init_b_face_num;    /* Initial border face numbering
                                    (if renumbering) */

  /* Periodictity features */

  cs_int_t  n_init_perio;         /* Number of initial periodicities
                                     (standard) */
  cs_int_t  n_transforms;         /* Number of transformations */

  fvm_periodicity_t  *periodicity;   /* parameters of each periodicity */

  /* Parallelism and/or periodic features */

  cs_mesh_halo_type_t  halo_type;  /* Halo type */

  cs_int_t   n_cells_with_ghosts;  /* Total number of cells on the local rank
                                      (n_cells + n_ghost_cells) */
  cs_int_t   n_ghost_cells;        /* Number of "ghost" cells */

  cs_mesh_halo_t  *halo;           /* Structure used to manage ghost cells */

  /* Extended neighborhood features */

  cs_int_t  *vtx_gcells_idx;   /* Index of the connectivity vertex -> cells
                                  Used to build the connectivity cell -> cells */

  cs_int_t  *vtx_gcells_lst;   /* Connectivity vertex -> cells. */

  cs_int_t  *cell_cells_idx;   /* "cell -> cells" connectivity index for
                                  extended halo. Only defined if extended
                                  neighborhood is built. */
  cs_int_t  *cell_cells_lst;   /* "cell -> cells" connectivity list for
                                  extended halo. Only defined if extended
                                  neighborhood is built. */

  /* Group and family features */

  cs_int_t    n_groups;    /* Number of groups */
  cs_int_t   *group_idx;   /* Starting index in the in group_lst */
  char       *group_lst;   /* List of group names */

  cs_int_t    n_families;          /* Number of families */
  cs_int_t    n_max_family_items;  /* Max. number of items for one family */
  cs_int_t   *family_item;         /* Family items */
  cs_int_t   *cell_family;         /* Cell family */
  cs_int_t   *b_face_family;       /* Border face family */

  fvm_group_class_set_t *class_defs;

  fvm_selector_t *select_cells;       /* Cell selection object */
  fvm_selector_t *select_i_faces;     /* Internal faces selection object */
  fvm_selector_t *select_b_faces;     /* Border faces selection object */

} cs_mesh_t ;

/* Structure used for building mesh structure */
/* ------------------------------------------ */

typedef struct {

  /* Periodic features */

  cs_int_t   *per_face_idx;    /* Index on periodicity for per_face_lst */

  cs_int_t   *per_face_lst;    /* Periodic faces list. For each couple,
                                  we have the local face number on local rank
                                  and the local face number on distant rank */

  cs_int_t   *per_rank_lst;    /* Remote ranks list. For each couple,
                                  we have the distant rank number. Exist
                                  only in case of parallelism. */

  /* Temporary features to define a periodicity. This kind of information
     is currently sent by the pre-processor. */

  cs_int_t                 perio_num;
  fvm_periodicity_type_t   perio_type;
  cs_real_t                translation[3];
  cs_real_t                invariant_point[3];
  cs_real_t                rotation_matrix[3][3];

} cs_mesh_builder_t ;

/*============================================================================
 * Static global variables
 *============================================================================*/

extern cs_mesh_t *cs_glob_mesh; /* Pointer to main mesh structure */

extern cs_mesh_builder_t  *cs_glob_mesh_builder; /* Pointer on builder mesh
                                                    structure */

/*============================================================================
 *  Public functions definition for API Fortran
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Recherche du numéro de groupe correspondant à un nom donné. Si le groupe
 * existe, le numéro renvoyé correspond au rang du groupe (commençant à 1)
 * dans la liste des groupes du maillage, et multiplié par -1. Cette
 * numérotation est celle utilisée dans le tableau IPRFML(NFML, NPRFML)
 * de description des familles.
 *
 * Si le groupe de nom indiqué n'existe pas, on renvoie 9999.
 *
 * Interface Fortran :
 *
 * FUNCTION NUMGRP (NOMGRP, LNGNOM)
 * ***************
 *
 * CHARACTER*       NOMGRP      : --> : Nom du groupe recherché
 * INTEGER          LNGNOM      : --> : Longueur du nom du groupe
 *----------------------------------------------------------------------------*/

cs_int_t CS_PROCF (numgrp, NUMGRP)
(
 const char       *const nomgrp,  /* --> Nom du groupe                        */
 const cs_int_t   *const lngnom   /* --> Longueur du nom                      */
 CS_ARGF_SUPP_CHAINE              /*     (arguments 'longueur' éventuels F77, */
                                  /*     inutilisés lors de l'appel mais      */
                                  /*     placés par de nombreux compilateurs) */
);


/*----------------------------------------------------------------------------
 * Sauvegarde des numérotations initiales des faces
 *
 * Remarque : lorsque des faces sont renumérotées, les connectivités et
 *            variables basées sur les faces sont toutes mises à jour, sauf
 *            les numéros globaux associés à ces faces (en cas de parallélisme)
 *            et conservés dans la structure maillage ('num_fac'et num_fbr') ;
 *            en effet, cette numérotation est destinée à se ramener à la
 *            numérotation initiale, notamment pour l'écriture et la lecture de
 *            fichiers suite indépendants du nombre de domaines en parallélisme
 *            ou d'une renumérotation vectorielle, et il est plus aisé de la
 *            manipuler si elle n'est pas modifiée.
 *
 * Interface Fortran :
 *
 * SUBROUTINE SAVNUM (IVECTV, IVECTB, INUMFI, INUMFB)
 * *****************
 *
 * INTEGER          IVECTV      : --> : Indicateur renum. faces internes
 * INTEGER          IVECTB      : --> : Indicateur renum. faces de bord
 * INTEGER          INUMFI      : --> : Table de renum. des faces internes
 * INTEGER          INUMFB      : --> : Table de renum. des faces de bord
 *----------------------------------------------------------------------------*/

void CS_PROCF (savnum, SAVNUM)
(
 const cs_int_t   *const ivectv,  /* --> Indicateur renum. faces internes     */
 const cs_int_t   *const ivectb,  /* --> Indicateur renum. faces de bord      */
 const cs_int_t   *const inumfi,  /* --> Table de renum. des faces internes   */
 const cs_int_t   *const inumfb   /* --> Table de renum. des faces de bord    */
);


/*=============================================================================
 * Public function definitions
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
 * Create an empty mesh builder structure
 *
 * returns:
 *   A pointer to a mesh builder structure.
 *----------------------------------------------------------------------------*/

cs_mesh_builder_t *
cs_mesh_builder_create(void);

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
 * Destroy a mesh builder structure
 *
 * mesh_builder     <->  pointer to a mesh structure
 *
 * returns:
 *  NULL pointer
 *----------------------------------------------------------------------------*/

cs_mesh_builder_t *
cs_mesh_builder_destroy(cs_mesh_builder_t  *mesh_builder);

/*----------------------------------------------------------------------------
 * Renumber vertices.
 *
 * We ensure:
 * If i < j then mesh->global_vtx_num[i] < mesh->global_vtx_num[j]
 * which is not insured by the initial numbering from the pre-processor.
 *
 * parameters:
 *   mesh      <->  pointer to mesh structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_order_vertices(cs_mesh_t  *const mesh);

/*----------------------------------------------------------------------------
 * Print mesh characteristics
 *
 * parameters:
 *   mesh         <-- pointer to mesh structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_info(const cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------
 * Compute global number of elements (cells, vertices, internal and border
 * faces) and sync cell family.
 *
 * parameters:
 *   mesh   <->  pointer to a cs_mesh_t structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_init_parall(cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------
 * Creation and initialization of halo structure.
 * Treatment of parallel and/or periodic halo for standard and extended
 * ghost cells according to halo type.
 *
 * parameters:
 *   mesh  <->  pointer to mesh structure.
 *----------------------------------------------------------------------------*/

void
cs_mesh_init_halo(cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------
 * Assign selectors to global mesh.
 *
 * Should be called once the mesh is fully built.
 *----------------------------------------------------------------------------*/

void
cs_mesh_init_selectors(void);

/*----------------------------------------------------------------------------
 * Dump of a mesh structure.
 *
 * parameters:
 *   mesh  <->  pointer to mesh structure.
 *----------------------------------------------------------------------------*/

void
cs_mesh_dump(const cs_mesh_t  *const mesh);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CS_MESH_H__ */
