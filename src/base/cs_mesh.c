/*============================================================================
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

/*============================================================================
 * Main structure associated to a mesh
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>
#include <bft_timer.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <fvm_order.h>
#include <fvm_parall.h>
#include <fvm_interface.h>
#include <fvm_periodicity.h>
#include <fvm_selector.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_halo.h"
#include "cs_mesh_halo.h"
#include "cs_perio.h"
#include "cs_mesh_quantities.h"
#include "cs_ext_neighborhood.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_mesh.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

#define CS_MESH_N_SUBS  5

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Pointer on the main mesh */

cs_mesh_t  *cs_glob_mesh = NULL;

/* Pointer on the temporary mesh used for building main mesh */

cs_mesh_builder_t  *cs_glob_mesh_builder = NULL;

/*============================================================================
 * Public FORTRAN function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Check necessity of extended mesh from FORTRAN options.
 *
 * Interface Fortran :
 *
 * SUBROUTINE HALTYP (IVOSET)
 * *****************
 *
 * INTEGER          IVOSET      : --> : Indicator of necessity of extended mesh
 *----------------------------------------------------------------------------*/

extern void
CS_PROCF (haltyp, HALTYP)(const cs_int_t   *ivoset);

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Update cell family array in case of parallelism.
 *
 * This function aims at copying main values from cells in halo (id between 1
 * and n_cells) to ghost cells on distant ranks (id between n_cells + 1 to
 * n_cells_with_halo).
 *
 * parameters:
 *   mesh  <->  pointer to mesh structure
 *----------------------------------------------------------------------------*/

static void
_sync_cell_fam(cs_mesh_t   *const mesh)
{
  cs_halo_t  *halo = mesh->halo;

  if (halo == NULL)
    return;

  bft_printf(_(" Synchronisation des familles des cellules\n"));

  cs_halo_sync_num(halo, CS_HALO_EXTENDED, mesh->cell_family);
}

/*----------------------------------------------------------------------------
 * Compute the minimum and the maximum of a vector (locally).
 *
 * parameters:
 *   n_vals    --> local number of elements
 *   var       --> pointer to vector
 *   min       <-- minimum
 *   max       <-- maximum
 *----------------------------------------------------------------------------*/

static void
_compute_local_minmax(cs_int_t        n_vals,
                      const cs_int_t  var[],
                      cs_int_t       *min,
                      cs_int_t       *max)
{
  cs_int_t  i;
  cs_int_t  _min = var[0], _max = var[0];

  for (i = 1; i < n_vals; i++) {
    _min = CS_MIN(_min, var[i]);
    _max = CS_MAX(_max, var[i]);
  }

  if (min != NULL)  *min = _min;
  if (max != NULL)  *max = _max;
}

/*----------------------------------------------------------------------------
 * Display the distribution of values of a integer vector.
 *
 * parameters:
 *   n_vals    --> local number of elements
 *   var       --> pointer to vector
 *----------------------------------------------------------------------------*/

static void
_display_histograms(cs_int_t        n_vals,
                    const cs_int_t  var[])
{
  cs_int_t  i, j, k, val_max, val_min;
  double step;

  cs_int_t count[CS_MESH_N_SUBS];
  cs_int_t n_steps = CS_MESH_N_SUBS;

  /* Compute local min and max */

  if (n_vals == 0) {
    bft_printf(_("    aucune valeur\n"));
    return;
  }

  val_max = var[0];
  val_min = var[0];
  _compute_local_minmax(n_vals, var, &val_min, &val_max);

  bft_printf(_("    valeur minimale =       %10d\n"), (int)val_min);
  bft_printf(_("    valeur maximale =       %10d\n\n"), (int)val_max);

  /* Define axis subdivisions */

  for (j = 0; j < n_steps; j++)
    count[j] = 0;

  if (val_max - val_min > 0) {

    if (val_max-val_min < n_steps)
      n_steps = CS_MAX(1, floor(val_max-val_min));

    step = (double)(val_max - val_min) / n_steps;

    /* Loop on values */

    for (i = 0; i < n_vals; i++) {

      /* Associated subdivision */

      for (j = 0, k = 1; k < n_steps; j++, k++) {
        if (var[i] < val_min + k*step)
          break;
      }
      count[j] += 1;

    }

    for (i = 0, j = 1; i < n_steps - 1; i++, j++)
      bft_printf("    %3d : [ %10d ; %10d [ = %10d\n",
                 i+1,
                 (int)(val_min + i*step),
                 (int)(val_min + j*step),
                 (int)(count[i]));

    bft_printf("    %3d : [ %10d ; %10d ] = %10d\n",
               n_steps,
               (int)(val_min + (n_steps - 1)*step),
               (int)val_max,
               (int)(count[n_steps - 1]));

  }

  else { /* if (val_max == val_min) */
    bft_printf("    %3d : [ %10d ; %10d ] = %10d\n",
               1, (int)(val_min), (int)val_max, (int)n_vals);
  }

}

/*----------------------------------------------------------------------------
 * Write a sum-up about halo features in listing
 *
 * parameters:
 *   mesh                   --> pointer to cs_mesh_t structure
 *   interface_time         --> time elapsed in interface build
 *   halo_time              --> time elapsed in halo build
 *   ext_neighborhood_time  --> time elapsed in extended neighborhood
 *                              connectivity building
 *----------------------------------------------------------------------------*/

static void
_print_halo_info(cs_mesh_t  *mesh,
                 double      interface_time,
                 double      halo_time,
                 double      ext_neighborhood_time)
{
  cs_halo_t  *halo = mesh->halo;

  cs_int_t  *rank_buffer = NULL;

  /* Sum-up of the computional times */

  bft_printf(_("\n Bilan des temps pour la création du halo\n\n"));

  if (mesh->n_domains > 1 || mesh->n_init_perio > 0)
    bft_printf(_("     Création de l'interface :                  %.3g s\n"),
               interface_time);

  if (mesh->halo_type == CS_HALO_EXTENDED)
    bft_printf(_("     Création de la connectivité étendue :      %.3g s\n"),
               ext_neighborhood_time);

  bft_printf(_("     Création du halo :                         %.3g s\n\n"),
             halo_time);


  bft_printf(_("     Temps total de création du halo :          %.3g s\n\n"),
             halo_time + interface_time + ext_neighborhood_time);
  bft_printf(" ----------------------------------------------------------\n\n");

  /* Sum-up ghost cell distribution */

  bft_printf(_(" Nombre de cellules standard :                         %d\n"),
             mesh->n_cells);

  if (mesh->n_domains > 1) {

    BFT_MALLOC(rank_buffer, mesh->n_domains, cs_int_t);

#if defined(_CS_HAVE_MPI)
    MPI_Allgather(&(mesh->n_cells), 1, CS_MPI_INT,
                  rank_buffer     , 1, CS_MPI_INT, cs_glob_base_mpi_comm);
#endif

    bft_printf(_("\n    Histogramme de la distribution du nombre de cellules"
                 " sur les rangs :\n\n"));

    _display_histograms(mesh->n_domains, rank_buffer);

  } /* End if n_domains > 1 */

  bft_printf("\n ----------------------------------------------------------\n");
  bft_printf_flush();

  bft_printf(_(" Nombre de cellules + cellules halo :                  %d\n\n"),
             mesh->n_cells_with_ghosts);

  if (mesh->n_domains > 1) {

#if defined(_CS_HAVE_MPI)
    MPI_Allgather(&(mesh->n_cells_with_ghosts), 1, CS_MPI_INT,
                  rank_buffer, 1, CS_MPI_INT, cs_glob_base_mpi_comm);
#endif

    bft_printf(_("\n    Histogramme de la distribution du nombre de cellules"
                 "\n    standard + halo (NCELET) sur les rangs :\n\n"));

    _display_histograms(mesh->n_domains, rank_buffer);

  } /* End if n_domains > 1 */

  bft_printf("\n ----------------------------------------------------------\n");
  bft_printf_flush();

  if (halo != NULL) {

    cs_int_t  n_std_ghost_cells = halo->n_elts[CS_HALO_STANDARD];

    bft_printf(_("\n Nombre de cellules fantômes locales :          %10d\n"),
               mesh->n_ghost_cells);
    bft_printf(_("     dans le voisinage standard:                %10d\n"),
               n_std_ghost_cells);
    bft_printf(_("     dans le voisinage etendu:                  %10d\n"),
               mesh->n_ghost_cells - n_std_ghost_cells);

    if (mesh->n_domains > 1) {

      cs_int_t  n_gcells = mesh->n_ghost_cells;

#if defined(_CS_HAVE_MPI)
      MPI_Allgather(&n_gcells  , 1, CS_MPI_INT,
                    rank_buffer, 1, CS_MPI_INT, cs_glob_base_mpi_comm);
#endif

      bft_printf(_("\n    Histogramme de la distribution des cellules"
                   " fantômes sur les rangs :\n\n"));

      _display_histograms(mesh->n_domains, rank_buffer);

    } /* If n_ranks > 1 */

    bft_printf("\n"
               " ----------------------------------------------------------\n");
    bft_printf_flush();

    /* Sum-up of the number of neighbors */

    bft_printf(_("\n Nombre de domaines voisins :         %d\n"),
               halo->n_c_domains);

    if (mesh->n_domains > 1) {

      cs_int_t  n_c_domains = halo->n_c_domains;

#if defined(_CS_HAVE_MPI)
      MPI_Allgather(&n_c_domains, 1, CS_MPI_INT,
                    rank_buffer , 1, CS_MPI_INT, cs_glob_base_mpi_comm);
#endif

      bft_printf(_("\n    Histogramme de la distribution du nombre de"
                   " voisins sur les rangs :\n\n"));

      _display_histograms(mesh->n_domains, rank_buffer);

    } /* If n_ranks > 1 */

    bft_printf("\n"
               " ----------------------------------------------------------\n");
    bft_printf_flush();

  } /* End if halo != NULL */

  if (mesh->n_domains > 1)
    BFT_FREE(rank_buffer);

}

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
)
{
  int  i, lngcmp ;

  char  *nomcmp = NULL;
  cs_mesh_t  *mesh = cs_glob_mesh;

  for (i = 0; i < mesh->n_groups; i++) {

    nomcmp = mesh->group_lst + (mesh->group_idx[i] - 1);
    lngcmp = strlen(nomcmp);

    if (lngcmp == *lngnom && strncmp(nomcmp, nomgrp, lngcmp) == 0)
      return (-(i + 1));

  }

  return -9999;

}

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
 * SUBROUTINE SAVNUM (IVECTI, IVECTB, INUMFI, INUMFB)
 * *****************
 *
 * INTEGER          IVECTI      : --> : Indicateur renum. faces internes
 * INTEGER          IVECTB      : --> : Indicateur renum. faces de bord
 * INTEGER          INUMFI      : --> : Table de renum. des faces internes
 * INTEGER          INUMFB      : --> : Table de renum. des faces de bord
 *----------------------------------------------------------------------------*/

void CS_PROCF (savnum, SAVNUM)
(
 const cs_int_t   *const ivecti,  /* --> Indicateur renum. faces internes     */
 const cs_int_t   *const ivectb,  /* --> Indicateur renum. faces de bord      */
 const cs_int_t   *const inumfi,  /* --> Table de renum. des faces internes   */
 const cs_int_t   *const inumfb   /* --> Table de renum. des faces de bord    */
)
{
  cs_int_t  i;

  cs_mesh_t  *mesh = cs_glob_mesh;

  /* Sauvegarde correspondance des faces internes -> faces internes initiales */

  if (*ivecti != 0) {

    BFT_MALLOC(mesh->init_i_face_num, mesh->n_i_faces, cs_int_t);

    for (i = 0; i < mesh->n_i_faces; i++)
      mesh->init_i_face_num[i] = inumfi[i];

  }

  /* Sauvegarde correspondance des faces de bord -> faces de bord initiales */

  if (*ivectb != 0) {

    BFT_MALLOC(mesh->init_b_face_num, mesh->n_b_faces, cs_int_t);

    for (i = 0; i < mesh->n_b_faces; i++)
      mesh->init_b_face_num[i] = inumfb[i];

  }

}

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
cs_mesh_create(void)
{
  cs_mesh_t  *mesh = NULL;

  BFT_MALLOC (mesh, 1, cs_mesh_t);

  /* General features */

  mesh->dim = 3;
  mesh->domain_num = cs_glob_base_rang + 1;
  mesh->n_domains = 0;

  /* Local dimensions */

  mesh->n_cells = 0;
  mesh->n_i_faces = 0;
  mesh->n_b_faces = 0;
  mesh->n_vertices = 0;
  mesh->i_face_vtx_connect_size = 0;
  mesh->b_face_vtx_connect_size = 0;

  /* Global dimensions */

  mesh->n_g_cells = 0;
  mesh->n_g_i_faces = 0;
  mesh->n_g_b_faces = 0;
  mesh->n_g_vertices = 0;

  /* Local structures */

  mesh->vtx_coord = NULL;
  mesh->i_face_cells = NULL;
  mesh->b_face_cells = NULL;
  mesh->i_face_vtx_idx = NULL;
  mesh->b_face_vtx_idx = NULL;
  mesh->i_face_vtx_lst = NULL;
  mesh->b_face_vtx_lst = NULL;

  /* Global numbering */

  mesh->global_cell_num = NULL;
  mesh->global_i_face_num = NULL;
  mesh->global_b_face_num = NULL;
  mesh->global_vtx_num = NULL;

  /* Initial face numbering (for exchange and restart) */

  mesh->init_i_face_num = NULL;
  mesh->init_b_face_num = NULL;

  /* Periodic features */

  mesh->n_init_perio = 0;
  mesh->n_transforms = 0;
  mesh->periodicity = NULL;

  /* Halo features */

  mesh->halo_type = CS_HALO_N_TYPES;
  mesh->n_ghost_cells = 0;
  mesh->n_cells_with_ghosts = 0;
  mesh->halo = NULL;

  mesh->vtx_gcells_idx = NULL;
  mesh->vtx_gcells_lst = NULL;
  mesh->cell_cells_idx = NULL;
  mesh->cell_cells_lst = NULL;

  /* Group and family features */

  mesh->n_groups = 0;
  mesh->group_idx = NULL;
  mesh->group_lst = NULL;

  mesh->n_max_family_items = 0;
  mesh->n_families = 0;

  mesh->family_item = NULL;
  mesh->cell_family = NULL;
  mesh->b_face_family = NULL;

  /* Selector features */

  mesh->class_defs = NULL;

  mesh->select_cells = NULL;
  mesh->select_i_faces = NULL;
  mesh->select_b_faces = NULL;

  return (mesh);
}

/*----------------------------------------------------------------------------
 * Create an empty mesh builder structure.
 *
 * returns:
 *   A pointer to a mesh builder structure
 *----------------------------------------------------------------------------*/

cs_mesh_builder_t *
cs_mesh_builder_create(void)
{
  cs_mesh_builder_t  *mesh_builder = NULL;

  BFT_MALLOC(mesh_builder, 1, cs_mesh_builder_t);

  mesh_builder->per_face_idx = NULL;
  mesh_builder->per_face_lst = NULL;
  mesh_builder->per_rank_lst = NULL;

  return mesh_builder;
}

/*----------------------------------------------------------------------------
 * Destroy a mesh structure.
 *
 * mesh       <->  pointer to a mesh structure
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

cs_mesh_t *
cs_mesh_destroy(cs_mesh_t  *mesh)
{
  BFT_FREE(mesh->vtx_coord);
  BFT_FREE(mesh->i_face_cells);
  BFT_FREE(mesh->b_face_cells);
  BFT_FREE(mesh->i_face_vtx_idx);
  BFT_FREE(mesh->b_face_vtx_idx);
  BFT_FREE(mesh->i_face_vtx_lst);
  BFT_FREE(mesh->b_face_vtx_lst);

  BFT_FREE(mesh->global_cell_num);
  BFT_FREE(mesh->global_i_face_num);
  BFT_FREE(mesh->global_b_face_num);
  BFT_FREE(mesh->global_vtx_num);

  BFT_FREE(mesh->init_i_face_num);
  BFT_FREE(mesh->init_b_face_num);

  BFT_FREE(mesh->group_idx);
  BFT_FREE(mesh->group_lst);

  BFT_FREE(mesh->family_item);
  BFT_FREE(mesh->cell_family);
  BFT_FREE(mesh->b_face_family);

  /* Free periodic structures */

  if (mesh->n_init_perio > 0)
    mesh->periodicity = fvm_periodicity_destroy(mesh->periodicity);

  if (mesh->halo_type == CS_HALO_EXTENDED) {
    BFT_FREE(mesh->vtx_gcells_idx);
    BFT_FREE(mesh->vtx_gcells_lst);
  }

  if (mesh->cell_cells_idx != NULL) {
    BFT_FREE(mesh->cell_cells_idx);
    BFT_FREE(mesh->cell_cells_lst);
  }

  /* Free halo structure */

  mesh->halo = cs_halo_destroy(mesh->halo);

  /* Free selection structures */

  if (mesh->n_groups > 0) {
    BFT_FREE(mesh->group_idx);
    BFT_FREE(mesh->group_lst);
  }

  if (mesh->select_cells != NULL)
    mesh->select_cells = fvm_selector_destroy(mesh->select_cells);
  if (mesh->select_i_faces != NULL)
    mesh->select_i_faces = fvm_selector_destroy(mesh->select_i_faces);
  if (mesh->select_b_faces != NULL)
    mesh->select_b_faces = fvm_selector_destroy(mesh->select_b_faces);

  /* Destroy group class set after selectors, who reference it */

  if (cs_glob_mesh->class_defs != NULL)
    cs_glob_mesh->class_defs
      = fvm_group_class_set_destroy(cs_glob_mesh->class_defs);

  BFT_FREE(mesh);

  return mesh;
}

/*----------------------------------------------------------------------------
 * Destroy a mesh builder structure
 *
 * mesh_builder     <->  pointer to a mesh structure
 *
 * returns:
 *  NULL pointer
 *----------------------------------------------------------------------------*/

cs_mesh_builder_t *
cs_mesh_builder_destroy(cs_mesh_builder_t  *mesh_builder)
{

  BFT_FREE(mesh_builder->per_face_idx);
  BFT_FREE(mesh_builder->per_face_lst);

  if (cs_glob_base_nbr > 1)
    BFT_FREE(mesh_builder->per_rank_lst);

  BFT_FREE(mesh_builder);

  return mesh_builder;
}

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
cs_mesh_order_vertices(cs_mesh_t  *const mesh)
{
  cs_int_t  i, j, size, dim, n_vertices;

  cs_int_t  *tmp_num = NULL;
  cs_real_t  *tmp_coord = NULL;
  fvm_lnum_t  *vertex_order = NULL;
  fvm_lnum_t  *vertex_renum = NULL;
  fvm_gnum_t  *g_vertex_num = NULL;

  assert(mesh != NULL);

  /* No treatment in serial mode */

  if (mesh->global_vtx_num == NULL)
    return;

  dim = mesh->dim;
  n_vertices = mesh->n_vertices;

  /* Compute the new vertex numbering */

  BFT_MALLOC(g_vertex_num, n_vertices, fvm_gnum_t);

  for (i = 0; i < n_vertices; i++)
    g_vertex_num[i] = mesh->global_vtx_num[i];

  vertex_order = fvm_order_local(NULL, g_vertex_num, (size_t)n_vertices);
  BFT_FREE(g_vertex_num);

  vertex_renum = fvm_order_local_renumbering(vertex_order, (size_t)n_vertices);
  BFT_FREE(vertex_order);

  /* Re-define face -> vertices connectivity with the new vertex numbering */

  if (mesh->n_i_faces > 0)
    for (i = 0, size = mesh->i_face_vtx_idx[mesh->n_i_faces] - 1; i < size; i++)
      mesh->i_face_vtx_lst[i] = vertex_renum[mesh->i_face_vtx_lst[i] - 1] + 1;

  if (mesh->n_b_faces > 0)
    for (i = 0, size = mesh->b_face_vtx_idx[mesh->n_b_faces] - 1; i < size; i++)
      mesh->b_face_vtx_lst[i] = vertex_renum[mesh->b_face_vtx_lst[i] - 1] + 1;

  /* Update coordinates */

  BFT_MALLOC(tmp_coord, n_vertices * dim, cs_real_t);

  for (i = 0; i < n_vertices; i++)
    for (j = 0; j < dim; j++)
      tmp_coord[vertex_renum[i]*dim + j] = mesh->vtx_coord[i*dim + j];

  memcpy(mesh->vtx_coord, tmp_coord, (n_vertices * dim * sizeof(cs_real_t)));
  BFT_FREE(tmp_coord);

  /* Update global numbering */

  BFT_MALLOC(tmp_num, n_vertices, cs_int_t);

  for (i = 0; i < n_vertices; i++)
    tmp_num[vertex_renum[i]] = mesh->global_vtx_num[i];

  memcpy(mesh->global_vtx_num, tmp_num, (n_vertices * sizeof(cs_int_t)));

  /* Free memory */

  BFT_FREE(tmp_num);
  BFT_FREE(vertex_renum);

}

/*----------------------------------------------------------------------------
 * Print mesh characteristics.
 *
 * parameters:
 *   mesh <-- pointer to mesh structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_info(const cs_mesh_t  *mesh)
{
  cs_int_t  i, vtx_id;

  cs_int_t  dim = mesh->dim;

  /* Allocation of the structures if this is not an update */

  if (mesh->n_g_vertices > 0) {

    cs_real_t  min_xyz[3] = { 1.e127,  1.e127,  1.e127};
    cs_real_t  max_xyz[3] = {-1.e127, -1.e127, -1.e127};

    for (vtx_id = 0 ; vtx_id < mesh->n_vertices ; vtx_id++) {

      for (i = 0 ; i < dim ; i++) {

        if (mesh->vtx_coord[vtx_id*dim + i] < min_xyz[i])
          min_xyz[i] = mesh->vtx_coord[vtx_id*dim + i];

        if (mesh->vtx_coord[vtx_id*dim + i] > max_xyz[i])
          max_xyz[i] = mesh->vtx_coord[vtx_id*dim + i];

      }

    } /* End of loop on vertices */


#if defined(_CS_HAVE_MPI)

    if (cs_glob_base_nbr > 1) {

      cs_real_t  g_min_xyz[3];
      cs_real_t  g_max_xyz[3];

      MPI_Allreduce(min_xyz, g_min_xyz, dim, CS_MPI_REAL, MPI_MIN,
                    cs_glob_base_mpi_comm);
      MPI_Allreduce(max_xyz, g_max_xyz, dim, CS_MPI_REAL, MPI_MAX,
                    cs_glob_base_mpi_comm);

      for (i = 0 ; i < dim ; i++) {
        min_xyz[i] = g_min_xyz[i];
        max_xyz[i] = g_max_xyz[i];
      }

    }

#endif

    bft_printf(_("\n"
                 " Coordonnées du maillage         minimale    et maximale\n"
                 "                       X : %14.7e %14.7e\n"
                 "                       Y : %14.7e %14.7e\n"
                 "                       Z : %14.7e %14.7e\n"),
               min_xyz[0], max_xyz[0], min_xyz[1], max_xyz[1],
               min_xyz[2], max_xyz[2]);

  }

}

/*----------------------------------------------------------------------------
 * Compute global number of elements (cells, vertices, internal and border
 * faces) and sync cell family.
 *
 * parameters:
 *   mesh  <->  pointer to a cs_mesh_t structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_init_parall(cs_mesh_t  *mesh)
{

#if defined(_CS_HAVE_MPI)

  cs_int_t  i;
  fvm_gnum_t  n_g_elts[4], max_elt_num[4];

  if (cs_glob_base_nbr <= 1)
    return;

  bft_printf(_("\n Définition globale du nombre d'éléments "
               "(cellules, sommets, faces...)\n"));

  /* Global dimensions of the mesh */

  max_elt_num[0] = mesh->n_cells;
  MPI_Allreduce(max_elt_num, n_g_elts, 1, FVM_MPI_GNUM, MPI_SUM,
                cs_glob_base_mpi_comm);

  max_elt_num[1] = 0;
  for (i = 0; i < mesh->n_i_faces; i++) {
    if (mesh->global_i_face_num[i] > max_elt_num[1])
      max_elt_num[1] = mesh->global_i_face_num[i];
  }

  max_elt_num[2] = 0;
  for (i = 0; i < mesh->n_b_faces; i++) {
    if (mesh->global_b_face_num[i] > max_elt_num[2])
      max_elt_num[2] = mesh->global_b_face_num[i];
  }

  max_elt_num[3] = 0;
  for (i = 0; i < mesh->n_vertices; i++) {
    if (mesh->global_vtx_num[i] > max_elt_num[3])
      max_elt_num[3] = mesh->global_vtx_num[i];
  }

  MPI_Allreduce(max_elt_num + 1, n_g_elts + 1, 3, FVM_MPI_GNUM, MPI_MAX,
                cs_glob_base_mpi_comm);

  mesh->n_g_cells = n_g_elts[0];
  mesh->n_g_i_faces = n_g_elts[1];
  mesh->n_g_b_faces = n_g_elts[2];
  mesh->n_g_vertices = n_g_elts[3];

#endif

  /* Sync cell family array (also in case of periodicity) */

  _sync_cell_fam(mesh);
}

/*----------------------------------------------------------------------------
 * Creation and initialization of halo structures.
 *
 * Treatment of parallel and/or periodic halos for standard and extended
 * ghost cells according to halo type requested by global options.
 *
 * parameters:
 *   mesh  <->  pointer to mesh structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_init_halo(cs_mesh_t  *mesh)
{
  int  n_periodic_lists, ivoset;
  fvm_lnum_t  i;
  double  t1, t2;
  double  halo_time = 0, interface_time = 0, ext_neighborhood_time = 0;

  int  *periodic_num = NULL;
  cs_int_t  *gcell_vtx_idx = NULL, *gcell_vtx_lst = NULL;
  fvm_lnum_t  *n_periodic_couples = NULL;
  fvm_gnum_t  *g_vertex_num = NULL;
  fvm_gnum_t  **periodic_couples = NULL;

  fvm_interface_set_t  *vertex_interfaces = NULL;

  const fvm_lnum_t  n_vertices = mesh->n_vertices;

  /* Choose the type of halo to build according to Fortran options.
     IMRGRA == 3 or 2 OR ITURB == 41 => CS_MESH_HALO_EXTENDED */

  CS_PROCF (haltyp, HALTYP) (&ivoset);

  /* Build halo */

  if (mesh->n_domains > 1 || mesh->n_init_perio > 0) {

    t1 = bft_timer_wtime();

    bft_printf("\n"
               " ----------------------------------------------------------\n");

    if (ivoset == 1) {

      bft_printf("\n Construction du halo avec voisinage étendu\n"
                 " ==========================================\n\n");

      mesh->halo_type = CS_HALO_EXTENDED;

      if (mesh->n_init_perio > 1) {

        bft_printf(_(" Composition des périodicités\n"));

        fvm_periodicity_combine(mesh->periodicity, 0);

      }

    }
    else {

      bft_printf("\n Construction du halo avec voisinage standard\n"
                 " ============================================\n\n");

      mesh->halo_type = CS_HALO_STANDARD;

    }

    /* Build purely parallel fvm_interface_set structure */

    if (mesh->n_domains == 1) {

      BFT_MALLOC(g_vertex_num, n_vertices, fvm_gnum_t);

      for (i = 0; i < n_vertices; i++)
        g_vertex_num[i] = (fvm_gnum_t)i+1;

    }
    else {

      assert(mesh->global_vtx_num != NULL);
      g_vertex_num = mesh->global_vtx_num;

    }

    if (mesh->n_init_perio > 0) {

      mesh->n_transforms = fvm_periodicity_get_n_transforms(mesh->periodicity);

      bft_printf(_(" Définition des couples périodiques\n"));
      bft_printf_flush();

      cs_perio_define_couples(&n_periodic_lists,
                              &periodic_num,
                              &n_periodic_couples,
                              &periodic_couples);

#if 0 /* For debugging purpose */
      for (i = 0; i < n_periodic_lists; i++) {

        cs_int_t  j;

        bft_printf("periodicity number: %d\n", periodic_num[i]);
        bft_printf("number of couples : %d\n", n_periodic_couples[i]);

        for (j = 0; j < n_periodic_couples[i]; j++)
          bft_printf("%4d --> %4d\n",
                     periodic_couples[i][2*j], periodic_couples[i][2*j + 1]);

      }

      fvm_periodicity_dump(mesh->periodicity);
#endif

    }

    bft_printf(_(" Création de l'interface\n"));
    bft_printf_flush();

    vertex_interfaces = fvm_interface_set_create(n_vertices,
                                                 NULL,
                                                 g_vertex_num,
                                                 mesh->periodicity,
                                                 n_periodic_lists,
                                                 periodic_num,
                                                 n_periodic_couples,
                       (const fvm_gnum_t **const)periodic_couples);

    if (mesh->n_domains == 1)
      BFT_FREE(g_vertex_num);
    BFT_FREE(periodic_num);
    BFT_FREE(n_periodic_couples);

    for (i = 0; i < mesh->n_init_perio; i++)
      BFT_FREE(periodic_couples[i]);
    BFT_FREE(periodic_couples);

#if 0 /* For debugging purposes */
    bft_printf("Dump de l'interface complete et finale\n");
    fvm_interface_set_dump(vertex_interfaces);
#endif

    t2 = bft_timer_wtime();
    interface_time = t2-t1;

    t1 = bft_timer_wtime();

    /* Creation of the cs_halo_t structure. */

    bft_printf(_(" Création de la structure halo\n"));
    bft_printf_flush();

    mesh->halo = cs_halo_create(vertex_interfaces);

    bft_printf(_(" Définition des halos\n"));
    bft_printf_flush();

    cs_mesh_halo_define(mesh,
                        vertex_interfaces,
                        &gcell_vtx_idx,
                        &gcell_vtx_lst);

    fvm_interface_set_destroy(vertex_interfaces);

    t2 = bft_timer_wtime();
    halo_time = t2-t1;

  } /* end if (mesh->n_domains > 1 || mesh->n_init_perio > 0) */

  /* Define a cell -> cells connectivity for the extended neighborhood
     if necessary */

  if (ivoset == 1) {

    t1 = bft_timer_wtime();
    bft_printf(_(" Définition des structures pour le voisinage étendu\n"));
    bft_printf_flush();

    cs_ext_neighborhood_define(mesh,
                               gcell_vtx_idx,
                               gcell_vtx_lst);

    bft_printf_flush();
    t2 = bft_timer_wtime();
    ext_neighborhood_time = t2-t1;

  }

  /* Free memory */

  BFT_FREE(gcell_vtx_idx);
  BFT_FREE(gcell_vtx_lst);

  /* Output for listing */

  if (mesh->halo_type != CS_HALO_N_TYPES)
    _print_halo_info(mesh,
                     interface_time,
                     halo_time,
                     ext_neighborhood_time);

  else if (ivoset == 1)
    bft_printf(_("\n Création de la connectivité étendue (%.3g s)\n"),
               ext_neighborhood_time);
}

/*----------------------------------------------------------------------------
 * Get the global number of ghost cells.
 *
 * parameters:
 *   mesh  -->  pointer to a mesh structure
 *
 * returns:
 *   global number of ghost cells
 *---------------------------------------------------------------------------*/

cs_int_t
cs_mesh_n_g_ghost_cells(cs_mesh_t  *mesh)
{
  cs_int_t  n_g_ghost_cells = 0;

  if (cs_glob_base_nbr == 1)
    n_g_ghost_cells = mesh->n_ghost_cells;

  else {

    assert(cs_glob_base_nbr > 1);

#if defined(_CS_HAVE_MPI)
    MPI_Allreduce(&(mesh->n_ghost_cells), &n_g_ghost_cells, 1, MPI_INT,
                  MPI_SUM, cs_glob_base_mpi_comm);
#endif

  }

  return n_g_ghost_cells;
}

/*----------------------------------------------------------------------------
 * Assign selectors to global mesh.
 *
 * Should be called once the mesh is fully built.
 *----------------------------------------------------------------------------*/

void
cs_mesh_init_selectors(void)
{
  int  i, j;
  int  grp_nbr, grp_num, grp_idx, color_nbr;

  int  *color = NULL;
  char **group = NULL;

  cs_glob_mesh->class_defs = fvm_group_class_set_create();

  /* Construction of the fvm_group_class structure */

  BFT_MALLOC(group, cs_glob_mesh->n_max_family_items, char*);
  BFT_MALLOC(color, cs_glob_mesh->n_max_family_items, int);

  for (i = 0; i < cs_glob_mesh->n_families; i++) {

    color_nbr = 0;
    grp_nbr  = 0;

    for (j = 0; j <  cs_glob_mesh->n_max_family_items; j++) {

      if (cs_glob_mesh->family_item[j * cs_glob_mesh->n_families + i] > 0){
        color[color_nbr++]
          = cs_glob_mesh->family_item[j *cs_glob_mesh->n_families + i];
      }
      else if (cs_glob_mesh->family_item[j * cs_glob_mesh->n_families + i] < 0) {
        /* Fortran formulation */
        grp_num = -cs_glob_mesh->family_item[j*cs_glob_mesh->n_families + i] -1;
        grp_idx = cs_glob_mesh->group_idx[grp_num];
        group[grp_nbr++] = cs_glob_mesh->group_lst + grp_idx -1;
      }

    }

    fvm_group_class_set_add(cs_glob_mesh->class_defs,
                            grp_nbr,
                            color_nbr,
                            (const char **)group,
                            color);

  } /* End of loop on families */

  BFT_FREE(group);
  BFT_FREE(color);

  /* Construction of the selectors */

  cs_glob_mesh->select_cells
    = fvm_selector_create(cs_glob_mesh->dim,
                          cs_glob_mesh->n_cells,
                          cs_glob_mesh->class_defs,
                          cs_glob_mesh->cell_family,
                          1,
                          cs_glob_mesh_quantities->cell_cen,
                          NULL);

  cs_glob_mesh->select_b_faces
    = fvm_selector_create(cs_glob_mesh->dim,
                          cs_glob_mesh->n_b_faces,
                          cs_glob_mesh->class_defs,
                          cs_glob_mesh->b_face_family,
                          1,
                          cs_glob_mesh_quantities->b_face_cog,
                          cs_glob_mesh_quantities->b_face_normal);

  cs_glob_mesh->select_i_faces
    = fvm_selector_create(cs_glob_mesh->dim,
                          cs_glob_mesh->n_i_faces,
                          cs_glob_mesh->class_defs,
                          NULL,
                          1,
                          cs_glob_mesh_quantities->i_face_cog,
                          cs_glob_mesh_quantities->i_face_normal);

}

/*----------------------------------------------------------------------------
 * Dump of a mesh structure.
 *
 * parameters:
 *   mesh  <->  pointer to mesh structure.
 *----------------------------------------------------------------------------*/

void
cs_mesh_dump(const cs_mesh_t  *const mesh)
{
  cs_int_t  i, j;

  bft_printf(_("\n\nDUMP OF THE MESH STRUCTURE: %p\n\n"),mesh);

  bft_printf(_("space dim :        %d\n"),mesh->dim);
  bft_printf(_("n_domains :        %d\n"),mesh->n_domains);
  bft_printf(_("domain_num:        %d\n"),mesh->domain_num);

  bft_printf(_("\nLocal dimensions:\n"));
  bft_printf(_("n_cells:                %d\n"),mesh->n_cells);
  bft_printf(_("n_cells_with_ghosts:        %d\n"),mesh->n_cells_with_ghosts);
  bft_printf(_("n_vertices:                %d\n"),mesh->n_vertices);
  bft_printf(_("n_i_faces:                %d\n"),mesh->n_i_faces);
  bft_printf(_("n_b_faces:                %d\n"),mesh->n_b_faces);

  bft_printf(_("\nGlobal dimensions:\n"));
  bft_printf(_("n_g_cells:        %d\n"),mesh->n_g_cells);
  bft_printf(_("n_g_vertices:        %d\n"),mesh->n_g_vertices);
  bft_printf(_("n_g_i_faces:        %d\n"),mesh->n_g_i_faces);
  bft_printf(_("n_g_b_faces:        %d\n"),mesh->n_g_b_faces);

  bft_printf(_("\n\n        --------"
               "        Vertices"
               "        --------\n\n"));

  bft_printf(_("\nVertex coordinates:\n"));
  for (i = 0; i < mesh->n_vertices; i++)
    bft_printf("<%3d > %10.3f        %10.3f        %10.3f\n",
               i+1, mesh->vtx_coord[3*i], mesh->vtx_coord[3*i+1],
               mesh->vtx_coord[3*i+2]);

  if (mesh->n_domains > 1) {
    bft_printf(_("\nGlobal vertex numbering:\n"));
    for (i = 0; i < mesh->n_vertices; i++)
      bft_printf("<%3d > %10u\n",
                 i+1, mesh->global_vtx_num[i]);
  }

  bft_printf(_("\n\n        ---------------------------"
               "        Internal faces connectivity"
               "        ---------------------------\n\n"));

  bft_printf(_("\nInternal faces -> Cells connectivity:\n"));
  for (i = 0; i < mesh->n_i_faces; i++)
    bft_printf("        < %d >         %5d    <----> %5d\n", i+1,
               mesh->i_face_cells[2*i], mesh->i_face_cells[2*i+1]);

  bft_printf(_("\nInternal faces -> vertices connectivity:\n"));
  for (i = 0; i < mesh->n_i_faces; i++) {
    bft_printf("        < %d >", i+1);
    for (j = mesh->i_face_vtx_idx[i]-1; j < mesh->i_face_vtx_idx[i+1]-1; j++)
      bft_printf("         %5d ",mesh->i_face_vtx_lst[j]);
    bft_printf("\n");
  }

  if (mesh->global_i_face_num != NULL) {

    bft_printf(_("\nInternal faces global numbering:\n"));
    for (i = 0; i < mesh->n_i_faces; i++)
      bft_printf("        < %d >         %12d",
                 i+1, mesh->global_i_face_num[i]);
    bft_printf("\n");

  }

  bft_printf(_("\n\n        -------------------------"
               "        Border faces connectivity"
               "        -------------------------\n\n"));

  bft_printf(_("\nBorder faces -> Cells connectivity:\n"));
  for (i = 0; i < mesh->n_b_faces; i++)
    bft_printf("        < %d >         %5d\n", i+1, mesh->b_face_cells[i]);

  bft_printf(_("\nBorder faces -> vertices connectivity:\n"));
  for (i = 0; i < mesh->n_b_faces; i++) {
    bft_printf("        < %d >", i+1);
    for (j = mesh->b_face_vtx_idx[i]-1; j < mesh->b_face_vtx_idx[i+1]-1; j++)
      bft_printf("         %5d ",mesh->b_face_vtx_lst[j]);
    bft_printf("\n");
  }

  bft_printf(_("\n\n        -------------------------"
               "        Cells"
               "        -------------------------\n\n"));

  if (mesh->global_cell_num != NULL) {

    bft_printf(_("\nCell global numbering:\n"));
    for (i = 0; i < mesh->n_cells; i++)
      bft_printf("        < %d >         %12d", i+1, mesh->global_cell_num[i]);
    bft_printf("\n");

  }

  bft_printf(_("\nNumber of families: %d\n"),mesh->n_families);
  bft_printf(_("Family of each cell:\n"));
  for (i = 0; i < mesh->n_cells_with_ghosts; i++)
    bft_printf("< %3d >         %5d\n", i+1, mesh->cell_family[i]);

  if (mesh->halo != NULL) {

    cs_halo_t  *halo = mesh->halo;

    bft_printf(_("\nHalo information: %p\n"), halo);

    bft_printf(_("n_c_domains:              %d\n"), halo->n_c_domains);
    bft_printf(_("n_ghost_cells:            %d\n"), mesh->n_ghost_cells);
    bft_printf(_("n_std_ghost_cells:        %d\n"), halo->n_elts[CS_HALO_STANDARD]);
    bft_printf(_("n_ext_ghost_cells:        %d\n"),
               halo->n_elts[CS_HALO_EXTENDED] - halo->n_elts[CS_HALO_STANDARD]);

    for (i = 0; i < halo->n_c_domains; i++) {

      bft_printf(_("\n\nRank id:        %d\n"
                   "Halo index start:        %d        end:        %d\n"
                   "Send index start:        %d        end:        %d\n"
                   "Send cell numbers:\n"),
                 halo->c_domain_rank[i],
                 halo->index[2*i], halo->index[2*i+2],
                 halo->send_index[2*i], halo->send_index[2*i+2]);
      for (j = halo->send_index[2*i]; j < halo->send_index[2*i+2]; j++)
        bft_printf("  %10d : %10d\n", j+1, halo->send_list[j]+1);

    } /* End of loop on the frontiers of halo */

    if (mesh->n_init_perio > 0 && halo->perio_lst != NULL) {

      const cs_int_t  n_c_domains = halo->n_c_domains;
      const cs_int_t  n_transforms = mesh->n_transforms;

      bft_printf(_("\n\nHalo's data in case of periodicity:\n"));
      bft_printf(_("n_transforms:                %d\n"),mesh->n_transforms);

      bft_printf("\nData in the standard halo\n");
      for (i = 0; i < n_transforms; i++)
        for (j = 0; j < n_c_domains; j++)
          bft_printf("< rank:%3d >< transform:%2d > start_idx: %5d"
                     "        n_elts: %5d\n",
                     halo->c_domain_rank[j], i,
                     halo->perio_lst[4*n_c_domains*i + 4*j],
                     halo->perio_lst[4*n_c_domains*i + 4*j+1]);

      bft_printf("\nData in the extended halo\n");
      for (i = 0; i < n_transforms; i++)
        for (j = 0; j < n_c_domains; j++)
          bft_printf("< rank:%3d >< transform:%2d >        "
                     "start_idx:  %5d, n_elts:  %5d\n",
                     halo->c_domain_rank[j], i,
                     halo->perio_lst[4*n_c_domains*i + 4*j+2],
                     halo->perio_lst[4*n_c_domains*i + 4*j+3]);

    } /* End if n_perio > 0 */

  } /* End if mesh->halo != NULL */

  if (mesh->cell_cells_idx != NULL) {

    bft_printf("\n\nCell -> cells connectivity for extended neighborhood\n\n");
    for (i = 0; i < mesh->n_cells; i++) {
      bft_printf("< cell n°:%3d>        ", i+1);
      for (j = mesh->cell_cells_idx[i]-1; j < mesh->cell_cells_idx[i+1]-1; j++)
        bft_printf("%d        ", mesh->cell_cells_lst[j]);
      bft_printf("\n");
    }

  }

  bft_printf(_("\n\nEND OF DUMP OF MESH STRUCTURE\n\n"));
  bft_printf_flush();

}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
