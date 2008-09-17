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
 * Fortran interface for reading data with "SolCom" specifications
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdarg.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_config.h>
#include <bft_mem.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <fvm_nodal.h>
#include <fvm_nodal_append.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_post.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_mesh_solcom.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Nombre de sommets, tétraèdres, pyramides, prismes, et hexaèdres */

static cs_int_t  cs_glob_nsom = 0;
static cs_int_t  cs_glob_ntetra = 0;
static cs_int_t  cs_glob_npyram = 0;
static cs_int_t  cs_glob_nprism = 0;
static cs_int_t  cs_glob_nhexae = 0;

/*============================================================================
 * Prototypes de fonctions privées
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Allocation de mémoire pour la lecture d'un fichier au format "SolCom" ;
 *----------------------------------------------------------------------------*/

static void cs_loc_maillage_solcom_alloc_mem
(
 cs_mesh_t             *const mesh,             /* <-> maillage associé       */
 cs_mesh_quantities_t  *const mesh_quantities   /* <-> grandeurs associés     */
);


/*----------------------------------------------------------------------------
 * Transfert d'une partie de la connectivité nodale à une structure FVM ;
 *----------------------------------------------------------------------------*/

static void cs_loc_maillage_solcom_ajoute
(
 fvm_nodal_t    *maillage_ext,            /* <-> maillage à compléter         */
 cs_int_t        nbr_elt,                 /* --> nombre d'élements à ajouter  */
 fvm_element_t   type,                    /* --> type de section à ajouter    */
 cs_int_t       *connect,                 /* --> connectivité à transférer    */
 cs_int_t       *cpt_elt_tot              /* <-> compteur total d'éléments    */
);


/*============================================================================
 * Prototypes de fonctions Fortran appellées
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Lecture des tableaux entités géométriques
 *----------------------------------------------------------------------------*/

void CS_PROCF (letgeo, LETGEO)
(
 const cs_int_t  *ndim,      /* --> dimension de l'espace                     */
 const cs_int_t  *ncelet,    /* --> nombre de cellules étendu                 */
 const cs_int_t  *ncel,      /* --> nombre de cellules                        */
 const cs_int_t  *nfac,      /* --> nombre de faces internes                  */
 const cs_int_t  *nfabor,    /* --> nombre de faces de bord                   */
 const cs_int_t  *nfml,      /* --> nombre de familles                        */
 const cs_int_t  *nprfml,    /* --> nombre de proprietes des familles         */
 const cs_int_t  *nnod,      /* --> nombre de sommets                         */
 const cs_int_t  *lndfac,    /* --> longueur de nodfac                        */
 const cs_int_t  *lndfbr,    /* --> longueur de nodfbr                        */
 const cs_int_t  *ntetra,    /* --> nombre de tetraèdres du maillage          */
 const cs_int_t  *npyram,    /* --> nombre de pyramides du maillage           */
 const cs_int_t  *nprism,    /* --> nombre de prismes du maillage             */
 const cs_int_t  *nhexae,    /* --> nombre d'hexaèdres du maillage            */
       cs_int_t  *inodal,    /* <-- indique si l'on doit lire la connectivite */
                             /*     nodale pour le post traitement            */
       cs_int_t   ifacel[],  /* <-- connect. faces internes / cellules        */
       cs_int_t   ifabor[],  /* <-- connect. faces de bord / cellules         */
       cs_int_t   ifmfbr[],  /* <-- liste des familles des faces bord         */
       cs_int_t   ifmcel[],  /* <-- liste des familles des cellules           */
       cs_int_t   iprfml[],  /* <-- liste des propriétés des familles         */
       cs_int_t   icotet[],  /* <-- connectivité nodale des tétraèdres        */
       cs_int_t   icopyr[],  /* <-- connectivité nodale des pyramides         */
       cs_int_t   icopri[],  /* <-- connectivité nodale des prismes           */
       cs_int_t   icohex[],  /* <-- connectivité nodale des hexaèdres         */
       cs_int_t   ipnfac[],  /* <-- rang ds nodfac 1er sommet faces int       */
       cs_int_t   nodfac[],  /* <-- numéro des sommets des faces int          */
       cs_int_t   ipnfbr[],  /* <-- rang ds nodfbr 1er sommt faces brd        */
       cs_int_t   nodfbr[],  /* <-- numéro des sommets des faces bord         */
       cs_real_t  xyzcen[],  /* <-- c.d.g. des cellules                       */
       cs_real_t  surfac[],  /* <-- surfaces des faces internes               */
       cs_real_t  surfbo[],  /* <-- surfaces des faces de bord                */
       cs_real_t  cdgfac[],  /* <-- c.d.g. des faces internes                 */
       cs_real_t  cdgfbo[],  /* <-- c.d.g. des faces de bord                  */
       cs_real_t  xyznod[]   /* <-- coordonnées des sommets                   */
);


/*============================================================================
 *  Fonctions publiques pour API Fortran
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Mise à jour des informations de dimensionnement maillage après la lecture
 * de l'entête du fichier de maillage en mode IFOENV = 0
 *
 * Interface Fortran :
 *
 * SUBROUTINE DIMGEO (NDIM  , NCELET, NCEL  , NFAC  , NFABOR, NSOM  ,
 * *****************
 *                    LNDFAC, LNDFBR, NFML  , NPRFML,
 *                    NTETRA, NPYRAM, NPRISM, NHEXAE )
 *
 * INTEGER          NDIM        : <-- : Dimension de l'espace (3)
 * INTEGER          NCELET      : <-- : Nombre d'éléments halo compris
 * INTEGER          NCEL        : <-- : Nombre d'éléments actifs
 * INTEGER          NFAC        : <-- : Nombre de faces internes
 * INTEGER          NFABOR      : <-- : Nombre de faces de bord
 * INTEGER          NSOM        : <-- : Nombre de sommets (optionnel)
 * INTEGER          LNDFAC      : <-- : Longueur de SOMFAC (optionnel)
 * INTEGER          LNDFBR      : <-- : Longueur de SOMFBR (optionnel)
 * INTEGER          NFML        : <-- : Nombre de familles des faces de bord
 * INTEGER          NPRFML      : <-- : Nombre de propriétés max par famille
 * INTEGER          NTETRA      : <-- : Nombre de tétraèdres
 * INTEGER          NPYRAM      : <-- : Nombre de pyramides
 * INTEGER          NPRISM      : <-- : Nombre de prismes
 * INTEGER          NHEXAE      : <-- : Nombre d'hexaèdres
 *----------------------------------------------------------------------------*/

void CS_PROCF (dimgeo, DIMGEO)
(
 const cs_int_t   *const ndim,    /* <-- dimension de l'espace                */
 const cs_int_t   *const ncelet,  /* <-- nombre d'éléments halo compris       */
 const cs_int_t   *const ncel,    /* <-- nombre d'éléments actifs             */
 const cs_int_t   *const nfac,    /* <-- nombre de faces internes             */
 const cs_int_t   *const nfabor,  /* <-- nombre de faces de bord              */
 const cs_int_t   *const nsom,    /* <-- nombre de sommets (optionnel)        */
 const cs_int_t   *const lndfac,  /* <-- longueur de somfac (optionnel)       */
 const cs_int_t   *const lndfbr,  /* <-- longueur de somfbr (optionnel)       */
 const cs_int_t   *const nfml,    /* <-- nombre de familles des faces de bord */
 const cs_int_t   *const nprfml,  /* <-- nombre de propriétés max par famille */
 const cs_int_t   *const ntetra,  /* <-- nombre de tétraèdres                 */
 const cs_int_t   *const npyram,  /* <-- nombre de pyramides                  */
 const cs_int_t   *const nprism,  /* <-- nombre de prismes                    */
 const cs_int_t   *const nhexae   /* <-- nombre d'hexaèdres                   */
)
{
  cs_mesh_t *mesh = cs_glob_mesh;

  mesh->dim = *ndim;

  mesh->n_cells = *ncel;
  mesh->n_i_faces = *nfac;
  mesh->n_b_faces = *nfabor;

  cs_glob_nsom = *nsom;

  if (*lndfac + *lndfbr > 0)
    mesh->n_vertices = *nsom;
  else
    mesh->n_vertices = 0;

  mesh->i_face_vtx_connect_size = *lndfac;
  mesh->b_face_vtx_connect_size = *lndfbr;

  mesh->n_cells_with_ghosts = *ncelet;

  assert (*ncelet == *ncel);

  mesh->n_g_cells = (fvm_gnum_t)mesh->n_cells;
  mesh->n_g_i_faces = (fvm_gnum_t)mesh->n_i_faces;
  mesh->n_g_b_faces = (fvm_gnum_t)mesh->n_b_faces;
  mesh->n_g_vertices = (fvm_gnum_t)mesh->n_vertices;

  mesh->n_max_family_items = *nprfml;
  mesh->n_families          = *nfml;

  cs_glob_ntetra = *ntetra;
  cs_glob_npyram = *npyram;
  cs_glob_nprism = *nprism;
  cs_glob_nhexae = *nhexae;

}


/*============================================================================
 * Fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Lecture d'un maillage au format "SolCom"
 *----------------------------------------------------------------------------*/

void cs_maillage_solcom_lit
(
 cs_mesh_t      *const mesh,      /* <-> maillage associé             */
 cs_mesh_quantities_t  *const mesh_quantities   /* <-> grandeurs associés           */
)
{

  cs_int_t   indic_nodal = 0;
  cs_int_t   cpt_elt_tot = 0;

  cs_real_t  *vtx_coord = NULL;
  cs_int_t   *connect_tetra = NULL;
  cs_int_t   *connect_pyram = NULL;
  cs_int_t   *connect_prism = NULL;
  cs_int_t   *connect_hexae = NULL;

  fvm_nodal_t  *maillage_ext = NULL;


  /* Allocations pour maillage principal */

  cs_loc_maillage_solcom_alloc_mem(mesh,
                                   mesh_quantities);


  /* Allocations pour maillage de post traitement
     quand on n'a pas de connectivité faces -> sommets */

  if (mesh->vtx_coord != NULL)
    vtx_coord = mesh->vtx_coord;

  else {

    BFT_MALLOC(vtx_coord, cs_glob_nsom * 3, cs_real_t);
    BFT_MALLOC(connect_tetra, cs_glob_ntetra * 4, cs_int_t);
    BFT_MALLOC(connect_pyram, cs_glob_npyram * 5, cs_int_t);
    BFT_MALLOC(connect_prism, cs_glob_nprism * 6, cs_int_t);
    BFT_MALLOC(connect_hexae, cs_glob_nhexae * 8, cs_int_t);

  }

  /* Lecture effective du corps du maillage */

  CS_PROCF (letgeo, LETGEO) (&(mesh->dim),
                             &(mesh->n_cells_with_ghosts),
                             &(mesh->n_cells),
                             &(mesh->n_i_faces),
                             &(mesh->n_b_faces),
                             &(mesh->n_families),
                             &(mesh->n_max_family_items),
                             &(cs_glob_nsom),
                             &(mesh->i_face_vtx_connect_size),
                             &(mesh->b_face_vtx_connect_size),
                             &cs_glob_ntetra,
                             &cs_glob_npyram,
                             &cs_glob_nprism,
                             &cs_glob_nhexae,
                             &indic_nodal,
                             mesh->i_face_cells,
                             mesh->b_face_cells,
                             mesh->b_face_family,
                             mesh->cell_family,
                             mesh->family_item,
                             connect_tetra,
                             connect_pyram,
                             connect_prism,
                             connect_hexae,
                             mesh->i_face_vtx_idx,
                             mesh->i_face_vtx_lst,
                             mesh->b_face_vtx_idx,
                             mesh->b_face_vtx_lst,
                             mesh_quantities->cell_cen,
                             mesh_quantities->i_face_normal,
                             mesh_quantities->b_face_normal,
                             mesh_quantities->i_face_cog,
                             mesh_quantities->b_face_cog,
                             vtx_coord);


  if (indic_nodal > 0) {

    /* Création directe du maillage de post-traitement
       lorsque l'on n'a pas de connectivité faces->sommets */

    maillage_ext = fvm_nodal_create(_("Fluid volume"), 3);

    if (cs_glob_ntetra > 0)
      cs_loc_maillage_solcom_ajoute(maillage_ext,
                                    cs_glob_ntetra,
                                    FVM_CELL_TETRA,
                                    connect_tetra,
                                    &cpt_elt_tot);

    if (cs_glob_npyram > 0)
      cs_loc_maillage_solcom_ajoute(maillage_ext,
                                    cs_glob_npyram,
                                    FVM_CELL_PYRAM,
                                    connect_pyram,
                                    &cpt_elt_tot);

    if (cs_glob_nprism > 0)
      cs_loc_maillage_solcom_ajoute(maillage_ext,
                                    cs_glob_nprism,
                                    FVM_CELL_PRISM,
                                    connect_prism,
                                    &cpt_elt_tot);

    if (cs_glob_nhexae > 0)
      cs_loc_maillage_solcom_ajoute(maillage_ext,
                                    cs_glob_nhexae,
                                    FVM_CELL_HEXA,
                                    connect_hexae,
                                    &cpt_elt_tot);

    fvm_nodal_transfer_vertices(maillage_ext,
                                vtx_coord);

    /* Transfert de la structure au post-traitement */

    cs_post_ajoute_maillage_existant(-1,
                                     maillage_ext,
                                     true);

  }
  else if (mesh->vtx_coord == NULL) {

    BFT_FREE(vtx_coord);
    BFT_FREE(connect_tetra);
    BFT_FREE(connect_pyram);
    BFT_FREE(connect_prism);
    BFT_FREE(connect_hexae);

  }

}


/*============================================================================
 * Fonctions privées
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Allocation de mémoire pour la lecture d'un fichier au format "SolCom" ;
 *----------------------------------------------------------------------------*/

static void cs_loc_maillage_solcom_alloc_mem
(
 cs_mesh_t             *const mesh,             /* <-> maillage associé       */
 cs_mesh_quantities_t  *const mesh_quantities   /* <-> grandeurs associés     */
)
{
  cs_int_t               nbr_elt;

  /* Allocation */
  /*------------*/

  /* Connectivités faces / cellules*/
  BFT_MALLOC(mesh->i_face_cells, mesh->n_i_faces * 2, cs_int_t);
  BFT_MALLOC(mesh->b_face_cells, mesh->n_b_faces, cs_int_t);

  /* CDG cellules (surdimensionnement usuel pour cellules fantômes,
     normalement inexistantes ici) */
  nbr_elt = mesh->dim * mesh->n_cells_with_ghosts;
  BFT_MALLOC(mesh_quantities->cell_cen, nbr_elt, cs_real_t);

  /* Surfaces faces */
  BFT_MALLOC(mesh_quantities->i_face_normal, mesh->dim * mesh->n_i_faces,
            cs_real_t);
  BFT_MALLOC(mesh_quantities->b_face_normal, mesh->dim * mesh->n_b_faces,
            cs_real_t);

  /* CDG faces */
  BFT_MALLOC(mesh_quantities->i_face_cog, mesh->dim * mesh->n_i_faces,
            cs_real_t);
  BFT_MALLOC(mesh_quantities->b_face_cog, mesh->dim * mesh->n_b_faces,
            cs_real_t);

  /* Familles faces de bord et cellules */
  BFT_MALLOC(mesh->b_face_family, mesh->n_b_faces, cs_int_t);
  BFT_MALLOC(mesh->cell_family, mesh->n_cells_with_ghosts, cs_int_t);

  /* Propriétés des familles */
  nbr_elt = mesh->n_families * mesh->n_max_family_items;
  BFT_MALLOC(mesh->family_item, nbr_elt, cs_int_t);


  if (mesh->n_vertices > 0) {

    /* Coordonnées des sommets */
    nbr_elt = mesh->dim * mesh->n_vertices;
    BFT_MALLOC(mesh->vtx_coord, nbr_elt, cs_real_t);

    /* Connectivités faces / sommets */
    BFT_MALLOC(mesh->i_face_vtx_idx, mesh->n_i_faces + 1, cs_int_t);
    BFT_MALLOC(mesh->i_face_vtx_lst, mesh->i_face_vtx_connect_size, cs_int_t);
    BFT_MALLOC(mesh->b_face_vtx_idx, mesh->n_b_faces + 1, cs_int_t);
    BFT_MALLOC(mesh->b_face_vtx_lst, mesh->b_face_vtx_connect_size, cs_int_t);

  }

}


/*----------------------------------------------------------------------------
 * Transfert d'une partie de la connectivité nodale à une structure FVM ;
 *----------------------------------------------------------------------------*/

static void cs_loc_maillage_solcom_ajoute
(
 fvm_nodal_t    *maillage_ext,            /* <-> maillage à compléter         */
 cs_int_t        nbr_elt,                 /* --> nombre d'élements à ajouter  */
 fvm_element_t   type,                    /* --> type de section à ajouter    */
 cs_int_t       *connect,                 /* --> connectivité à transférer    */
 cs_int_t       *cpt_elt_tot              /* <-> compteur total d'éléments    */
)
{

  cs_int_t   ind_elt = 0;
  cs_int_t   cpt_elt = 0;

  fvm_lnum_t   *parent_element_num = NULL;

  if (nbr_elt == 0)
    return;

  /* Numéros d'éléments parents à si nécessaire */

  if (cpt_elt > 0) {
    BFT_MALLOC(parent_element_num, nbr_elt, fvm_lnum_t);
    for (ind_elt = 0 ; ind_elt < nbr_elt ; ind_elt++)
      parent_element_num[ind_elt] = ind_elt + (*cpt_elt_tot + 1);
  }


  /* Transfert de la connectivité et des numéros de parents */

  fvm_nodal_append_by_transfer(maillage_ext,
                               nbr_elt,
                               type,
                               NULL,
                               NULL,
                               NULL,
                               connect,
                               parent_element_num);

  *cpt_elt_tot = *cpt_elt_tot + nbr_elt;

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
