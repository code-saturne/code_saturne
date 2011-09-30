/*============================================================================
 *  Définition de la fonction
 *   de lecture d'un fichier de maillage au format Comet
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

/*============================================================================
 *                                 Visibilité
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie standard C
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/

#include "ecs_def.h"
#include "ecs_elt_typ_liste.h"
#include "ecs_file.h"
#include "ecs_mem.h"
#include "ecs_tab.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles des paquetages visibles
 *----------------------------------------------------------------------------*/

#include "ecs_descr.h"
#include "ecs_descr_chaine.h"
#include "ecs_maillage.h"
#include "ecs_maillage_priv.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_maillage_pre.h"


/*----------------------------------------------------------------------------
 *  Fichier  `include' du  paquetage courant associe au fichier courant
 *----------------------------------------------------------------------------*/

#include "ecs_pre_comet.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' privés   du  paquetage courant
 *----------------------------------------------------------------------------*/


/*============================================================================
 *                    Déclaration de paramètres et macros
 *============================================================================*/


/*============================================================================
 *                  Définition de structures locales
 *============================================================================*/

/* Définition de noeuds */

typedef struct {
  ecs_int_t        nbr_noeuds;  /* Nombre de noeuds */
  ecs_int_t       *id_noeud;    /* Labels des noeuds */
  ecs_coord_t     *coord;       /* Coordonnées (entrelacées) */
} ecs_loc_noeuds_comet_t;


/* Définition de faces (bord et interne) */

typedef struct {
  ecs_int_t     nbr_fac;        /* Nombre de faces */
  ecs_int_t    *nbr_n;          /* Nbr. sommets/faces  */
  ecs_int_t     taille_connect; /* taille de la connectivite  */
  ecs_int_t    *connect;        /* Connectivité sommets */
  ecs_int_t    *icel1;          /* Cellule sur le cote positif de la face */
  ecs_int_t    *icel2;          /* Cellule sur le cote négatif de la face */
  ecs_int_t    *icoul;          /* Couleur des faces */
} ecs_loc_faces_comet_t;


/* Définition d'elements */

typedef struct {
  ecs_int_t        nbr_cel;    /* Nombre de faces */
  ecs_int_t       *id;         /* Couleur des faces */
  ecs_int_t       *icoul;      /* couleur des élements */
} ecs_loc_cels_comet_t;


typedef enum {
  VERTEX = 1,
  SHELL,
  BOUNDARY,
  INTERNAL_FACE,
  CELL,
  VERTEX_MAP,
  SHELL_MAP,
  BOUNDARY_MAP,
  INTERNAL_FACE_MAP,
  CELL_MAP,
  VERTEX_DOUBLE,
  INTERFACE,
  SCALAR
} ecs_loc_comet_typ_sec_t;


/*============================================================================
 *                              Fonctions privées
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Lecture d'une chaîne de caractères dans un fichier xdr
 *----------------------------------------------------------------------------*/

static char *
ecs_pre_comet__lit_chaine_xdr(ecs_file_t  *fic)
{
  ecs_int_t l_chaine;
  ecs_int_t nb_rec;
  ecs_int_t l_align;

  char * chaine_lue;
  char * chaine;
  int32_t *type;


  /*Xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


  /* Longueur de la chaîne */

  ECS_MALLOC(type, 1, int32_t);
  ecs_file_read(type, 4, 1, fic);

  l_chaine = (ecs_int_t) *type;

  ECS_FREE(type);

  ECS_MALLOC(chaine_lue, l_chaine + 1, char);

  ecs_file_read(chaine_lue, 1, l_chaine, fic);

  chaine_lue[l_chaine] = '\0';


  /* Alignement de 4 octets */

  for (l_align = l_chaine; l_align % 4 != 0; l_align++);

  nb_rec = l_align - l_chaine;

  if (nb_rec != 0) {
    ECS_MALLOC(chaine, nb_rec, char);
    ecs_file_read(chaine, 1, nb_rec, fic);
    ECS_FREE(chaine);
  }

  return chaine_lue;
}

/*----------------------------------------------------------------------------
 *  Lecture d'une section vertex et concaténation avec les éventuels sommets
 *   appartenant aux sections vertex préalablement lues
 *----------------------------------------------------------------------------*/

static void
ecs_pre_comet__lit_vertex(ecs_file_t               *fic,
                          ecs_loc_noeuds_comet_t  **noeuds,
                          ecs_int_t                *taille_noeuds,
                          bool                      is_little_endian)
{
  ecs_int_t nbr_noeuds_loc;
  ecs_int_t nbr_noeuds_ini;
  ecs_int_t nbr_noeuds_total;
  ecs_coord_t *coord_loc;
  ecs_int_t *id_loc;

  ecs_int_t taille_rec;
  char *tab_lec;
  ecs_int_t inoeud;

  int32_t *type = NULL;

  int couleur;

  /*Xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


  /* Lecture */

  ECS_MALLOC(type, 2, int32_t);

  ecs_file_read(type, 4, 2, fic);
  nbr_noeuds_loc = (ecs_int_t ) *type;
  couleur = *(type + 1);

  printf(_("\n    Reading a vertices section\n"));
  printf(_("         %10d vertices "), (int)nbr_noeuds_loc);
  printf(_("of color %4d\n\n"), couleur);

  ecs_file_read(type, 4, 1, fic);

  ECS_FREE(type);

  /* allocation */

  if (*taille_noeuds == 0) {
    *taille_noeuds = 1;
    ECS_REALLOC(*noeuds, *taille_noeuds, ecs_loc_noeuds_comet_t);
    nbr_noeuds_ini = 0;
    (*noeuds)->nbr_noeuds = nbr_noeuds_loc;
    ECS_MALLOC((*noeuds)->coord, 3 * nbr_noeuds_loc, ecs_coord_t);
    ECS_MALLOC((*noeuds)->id_noeud, nbr_noeuds_loc, ecs_int_t);
  }

  else {
    nbr_noeuds_ini = (*noeuds)->nbr_noeuds;
    nbr_noeuds_total = nbr_noeuds_ini + nbr_noeuds_loc;
    (*noeuds)->nbr_noeuds = nbr_noeuds_total;
    ECS_REALLOC((*noeuds)->coord, 3 * nbr_noeuds_total, ecs_coord_t);
    ECS_REALLOC((*noeuds)->id_noeud, nbr_noeuds_total, ecs_int_t);
  }

  coord_loc = (*noeuds)->coord + 3 * nbr_noeuds_ini;
  id_loc = (*noeuds)->id_noeud + nbr_noeuds_ini;


  /* Lecture de la section */

  taille_rec = 4 + (4 * 3); /* 1 entier + 3 réels */

  ECS_MALLOC(tab_lec, nbr_noeuds_loc * taille_rec, char);

  ecs_file_read(tab_lec, 1, nbr_noeuds_loc * taille_rec, fic);


  /* Test si le système est little/big endian */

  if (is_little_endian == true)
    ecs_file_swap_endian(tab_lec, tab_lec, 4, 4 * nbr_noeuds_loc);


  /* Coordonnées + id */

  for (inoeud = 0; inoeud < nbr_noeuds_loc; inoeud ++) {

    id_loc[inoeud]
      = (ecs_int_t) *((int32_t *)(tab_lec + taille_rec*inoeud));

    coord_loc[inoeud*3]
      = (ecs_coord_t) *((float *) (tab_lec + taille_rec*inoeud + 4));
    coord_loc[inoeud*3 + 1]
      = (ecs_coord_t) *((float *) (tab_lec + taille_rec*inoeud + 8));
    coord_loc[inoeud*3 + 2]
      = (ecs_coord_t) *((float *) (tab_lec + taille_rec*inoeud + 12));

  }

  ECS_FREE(tab_lec);

}

/*----------------------------------------------------------------------------
 *  Lecture d'une section vertex dont les coordonnées sont en double précision
 *   et concaténation avec les éventuels sommets
 *   appartenant aux sections vertex préalablement lues
 *----------------------------------------------------------------------------*/

static void
ecs_pre_comet__lit_vertex_dbl(ecs_file_t               *fic,
                              ecs_loc_noeuds_comet_t  **noeuds,
                              ecs_int_t                *taille_noeuds,
                              bool                      is_little_endian)
{
  ecs_int_t nbr_noeuds_loc;
  ecs_int_t nbr_noeuds_ini;
  ecs_int_t nbr_noeuds_total;
  ecs_coord_t *coord_loc;
  ecs_int_t *id_loc;

  ecs_int_t taille_rec;
  char *tab_lec;
  ecs_int_t inoeud;

  int32_t *type = NULL;

  int couleur;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Lecture */

  ECS_MALLOC(type, 2, int32_t);

  ecs_file_read(type, 4, 2, fic);
  nbr_noeuds_loc = (ecs_int_t ) *type;
  couleur = *(type + 1);

  printf(_("\n    Reading vertices with double precision\n"));

  printf(_("         %10d vertices "), (int)nbr_noeuds_loc);
  printf(_("of color %4d\n\n"), couleur);

  ecs_file_read(type, 4, 1, fic);

  ECS_FREE(type);


  /* allocation */

  if (*taille_noeuds == 0) {

    *taille_noeuds = 1;
    ECS_REALLOC(*noeuds, *taille_noeuds, ecs_loc_noeuds_comet_t);
    nbr_noeuds_ini = 0;
    (*noeuds)->nbr_noeuds = nbr_noeuds_loc;
    ECS_MALLOC((*noeuds)->coord, 3 * nbr_noeuds_loc, ecs_coord_t);
    ECS_MALLOC((*noeuds)->id_noeud, nbr_noeuds_loc, ecs_int_t);
  }

  else {

    nbr_noeuds_ini = (*noeuds)->nbr_noeuds;
    nbr_noeuds_total = nbr_noeuds_ini + nbr_noeuds_loc;
    (*noeuds)->nbr_noeuds = nbr_noeuds_total;
    ECS_REALLOC((*noeuds)->coord, 3 * nbr_noeuds_total, ecs_coord_t);
    ECS_REALLOC((*noeuds)->id_noeud, nbr_noeuds_total, ecs_int_t);
  }

  coord_loc = (*noeuds)->coord + 3 * nbr_noeuds_ini;
  id_loc = (*noeuds)->id_noeud + nbr_noeuds_ini;


  /* Lecture de la section */

  taille_rec = 4 + (8 * 3); /* 1 entier + 3 réels */

  ECS_MALLOC(tab_lec, nbr_noeuds_loc * taille_rec, char);

  ecs_file_read(tab_lec,
                1,
                nbr_noeuds_loc * taille_rec,
                fic);


  /* Coordonnées + id */

  for (inoeud = 0; inoeud < nbr_noeuds_loc; inoeud ++) {
   if (is_little_endian == true)
     ecs_file_swap_endian(tab_lec + taille_rec*inoeud,
                          tab_lec + taille_rec*inoeud,
                          4,
                          1);
   id_loc[inoeud]
     = (ecs_int_t) *((int32_t *)(tab_lec + taille_rec*inoeud));

   if (is_little_endian == true)
     ecs_file_swap_endian(tab_lec + taille_rec*inoeud + 4,
                          tab_lec + taille_rec*inoeud + 4,
                          8,
                          3);


    /* Vérification que ecs_coord_t est bien un double */

   assert (sizeof(ecs_coord_t) == 8 );

   coord_loc[inoeud*3    ]
     = *((ecs_coord_t *)(tab_lec + taille_rec*inoeud + 4));
   coord_loc[inoeud*3 + 1]
     = *((ecs_coord_t *)(tab_lec + taille_rec*inoeud + 12));
   coord_loc[inoeud*3 + 2]
     = *((ecs_coord_t *)(tab_lec + taille_rec*inoeud + 20));

   ECS_FREE(tab_lec);
  }

}

/*----------------------------------------------------------------------------
 *  Lecture d'une section internal face et concaténation avec les éventuels
 *  éléments appartenant aux sections internal face préalablement lues
 *----------------------------------------------------------------------------*/

static void
ecs_pre_comet__lit_face(ecs_file_t              *fic,
                        ecs_loc_faces_comet_t  **faces,
                        ecs_int_t               *taille_faces,
                        int                      version)
{
  ecs_int_t nbr_fac_loc;
  ecs_int_t nbr_fac_ini;
  ecs_int_t nbr_fac_total;
  ecs_int_t *icoul_loc;
  ecs_int_t *icel1_loc;
  ecs_int_t *icel2_loc;
  ecs_int_t *nbr_n_loc;
  ecs_int_t taille_connect_loc;

  ecs_int_t ifac;
  ecs_int_t isom;
  ecs_int_t cpt_som;
  ecs_int_t max_som;
  ecs_int_t max_som_fac;
  ecs_int_t cpt_som_total;

  int couleur;

  int32_t *tab_lec = NULL;
  int32_t *int_lec = NULL;


  /*Xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


  /* Lecture */

  ECS_MALLOC(int_lec, 3, int32_t);

  ecs_file_read(int_lec, 4, 3, fic);
  nbr_fac_loc = (ecs_int_t ) *int_lec;
  couleur = *(int_lec + 1);

  printf(_("\n    Reading an internal faces section\n"));
  printf(_("         %10d faces "), (int)nbr_fac_loc);
  printf(_("of color %4d\n\n"), couleur);

  /* allocation et approximation du nombre de sommets par face
     pour diminuer le nombre de réallocations */

  if (*taille_faces == 0) {
    *taille_faces = 1;
    ECS_REALLOC(*faces, *taille_faces, ecs_loc_faces_comet_t);
    nbr_fac_ini = 0;
    (*faces)->nbr_fac = nbr_fac_loc;
    (*faces)->taille_connect = 0;
    max_som = 4 * nbr_fac_loc;
    ECS_MALLOC((*faces)->icoul, nbr_fac_loc, ecs_int_t);
    ECS_MALLOC((*faces)->icel1, nbr_fac_loc, ecs_int_t);
    ECS_MALLOC((*faces)->icel2, nbr_fac_loc, ecs_int_t);
    ECS_MALLOC((*faces)->nbr_n, nbr_fac_loc, ecs_int_t);
    ECS_MALLOC((*faces)->connect, max_som, ecs_int_t);
  }

  else {
    nbr_fac_ini = (*faces)->nbr_fac;
    nbr_fac_total = nbr_fac_ini + nbr_fac_loc;
    (*faces)->nbr_fac = nbr_fac_total;
    max_som = 4 * nbr_fac_loc + (*faces)->taille_connect;
    ECS_REALLOC((*faces)->icoul, nbr_fac_total, ecs_int_t);
    ECS_REALLOC((*faces)->icel1, nbr_fac_total, ecs_int_t);
    ECS_REALLOC((*faces)->icel2, nbr_fac_total, ecs_int_t);
    ECS_REALLOC((*faces)->nbr_n, nbr_fac_total, ecs_int_t);
    ECS_REALLOC((*faces)->connect, max_som, ecs_int_t);
  }

  icoul_loc = (*faces)->icoul;
  icel1_loc = (*faces)->icel1;
  icel2_loc = (*faces)->icel2;
  nbr_n_loc = (*faces)->nbr_n;
  taille_connect_loc = (*faces)->taille_connect;

  ECS_REALLOC(int_lec, 6 , int32_t);

  cpt_som = 0;


  /* Allocation du tableau temporaire */

  max_som_fac = 8;
  ECS_MALLOC(tab_lec, max_som_fac , int32_t);


  /* lecture face par face */

  for (ifac = 0; ifac < nbr_fac_loc; ifac ++) {


    if (version <= 1) {

      ecs_file_read(int_lec, 4, 4, fic);

      /* icel1 */

      icel1_loc[nbr_fac_ini + ifac] = (ecs_int_t) *(int_lec + 1);

      /* icel2 */

      icel2_loc[nbr_fac_ini + ifac] = (ecs_int_t) *(int_lec + 2);

      /* nombres de sommets de la face */

      nbr_n_loc[nbr_fac_ini + ifac] = (ecs_int_t) *(int_lec + 3);

    }

    else {

      ecs_file_read(int_lec, 4, 6, fic);

      /* icel1 */

      icel1_loc[nbr_fac_ini + ifac] = (ecs_int_t) *(int_lec + 1);

      /* icel2 */

      icel2_loc[nbr_fac_ini + ifac] = (ecs_int_t) *(int_lec + 3);

      /* nombres de sommets de la face */

      nbr_n_loc[nbr_fac_ini + ifac] = (ecs_int_t) *(int_lec + 5);

    }

    icoul_loc[nbr_fac_ini + ifac] = couleur;

    /* connectivité */

    cpt_som_total =   taille_connect_loc + cpt_som
                    + nbr_n_loc[nbr_fac_ini + ifac];

    if (max_som < cpt_som_total) {
      ECS_REALLOC((*faces)->connect, 2 * max_som, ecs_int_t);
      max_som = 2 * max_som;
    }

    if (max_som_fac < nbr_n_loc[nbr_fac_ini + ifac]) {
      ECS_REALLOC(tab_lec, max_som_fac * 2, int32_t);
      max_som_fac = max_som_fac * 2;
    }

    ecs_file_read(tab_lec, 4, nbr_n_loc[nbr_fac_ini + ifac], fic);

    for (isom = 0; isom < nbr_n_loc[nbr_fac_ini + ifac]; isom ++)
      (*faces)->connect[taille_connect_loc + cpt_som + isom]
        = (ecs_int_t) tab_lec[isom];


    /* mise à jour du compteur */

    cpt_som += nbr_n_loc[nbr_fac_ini + ifac];

  }


  /* Mise à jour de la taille de la connectivité */

  (*faces)->taille_connect = taille_connect_loc + cpt_som;
  ECS_REALLOC((*faces)->connect,  taille_connect_loc + cpt_som, ecs_int_t);

  ECS_FREE(int_lec);
  ECS_FREE(tab_lec);

}

/*----------------------------------------------------------------------------
 *  Lecture d'une section boundary et concatenation avec les eventuelles
 *  faces de bord appartenant aux sections boundary prealablement lues
 *----------------------------------------------------------------------------*/

static void
ecs_pre_comet__lit_boundary(ecs_file_t              *fic,
                            ecs_loc_faces_comet_t  **faces,
                            ecs_int_t               *taille_faces,
                            int                      version)
{
  ecs_int_t nbr_fac_loc;
  ecs_int_t nbr_fac_ini;
  ecs_int_t nbr_fac_total;
  ecs_int_t *icoul_loc;
  ecs_int_t *icel1_loc;
  ecs_int_t *icel2_loc;
  ecs_int_t *nbr_n_loc;
  ecs_int_t taille_connect_loc;

  ecs_int_t ifac;
  ecs_int_t isom;
  ecs_int_t cpt_som;
  ecs_int_t max_som;
  ecs_int_t max_som_fac;
  ecs_int_t cpt_som_total;

  int couleur;

  int32_t *tab_lec = NULL;
  int32_t *int_lec = NULL;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Lecture */

  ECS_MALLOC(int_lec, 3, int32_t);

  ecs_file_read(int_lec, 4, 3, fic);
  nbr_fac_loc = (ecs_int_t ) *int_lec;
  couleur = *(int_lec + 1);

  printf(_("\n    Reading a boundary faces section\n"));
  printf(_("         %10d boundary faces "), (int)nbr_fac_loc);
  printf(_("of color %4d\n\n"), couleur);

  /* allocation */

  if (*taille_faces == 0) {
    *taille_faces = 1;
    ECS_REALLOC(*faces, *taille_faces, ecs_loc_faces_comet_t);
    nbr_fac_ini = 0;
    (*faces)->nbr_fac = nbr_fac_loc;
    (*faces)->taille_connect = 0;
    ECS_MALLOC((*faces)->icoul, nbr_fac_loc, ecs_int_t);
    ECS_MALLOC((*faces)->icel1, nbr_fac_loc, ecs_int_t);
    ECS_MALLOC((*faces)->icel2, nbr_fac_loc, ecs_int_t);
    ECS_MALLOC((*faces)->nbr_n, nbr_fac_loc, ecs_int_t);
    max_som = 4 * nbr_fac_loc;
    ECS_MALLOC((*faces)->connect, max_som, ecs_int_t);
  }

  else {
    nbr_fac_ini = (*faces)->nbr_fac;
    nbr_fac_total = nbr_fac_ini + nbr_fac_loc;
    (*faces)->nbr_fac = nbr_fac_total;
    max_som = 4 * nbr_fac_loc + (*faces)->taille_connect;
    ECS_REALLOC((*faces)->icoul, nbr_fac_total, ecs_int_t);
    ECS_REALLOC((*faces)->icel1, nbr_fac_total, ecs_int_t);
    ECS_REALLOC((*faces)->icel2, nbr_fac_total, ecs_int_t);
    ECS_REALLOC((*faces)->nbr_n, nbr_fac_total, ecs_int_t);
    ECS_REALLOC((*faces)->connect, max_som, ecs_int_t);
  }

  icoul_loc = (*faces)->icoul;
  icel1_loc = (*faces)->icel1;
  icel2_loc = (*faces)->icel2;
  nbr_n_loc = (*faces)->nbr_n;
  taille_connect_loc = (*faces)->taille_connect;

  ECS_REALLOC(int_lec, 6 , int32_t);

  /* lecture face par face */

  cpt_som = 0;

  /* Allocation du tableau temporaire */

  max_som_fac = 8;
  ECS_MALLOC(tab_lec, max_som_fac , int32_t);

  for (ifac = 0; ifac < nbr_fac_loc; ifac ++) {

    if (version <= 1) {

      ecs_file_read(int_lec, 4, 3, fic);

      /* icel1 */

      icel1_loc[nbr_fac_ini + ifac] = -1;

      /* icel2 */

      icel2_loc[nbr_fac_ini + ifac] = (ecs_int_t) *(int_lec + 1);

      /* nombres de sommets de la face */

      nbr_n_loc[nbr_fac_ini + ifac] = (ecs_int_t) *(int_lec + 2);

    }

    else {

      ecs_file_read(int_lec, 4, 6, fic);

      /* icel1 */

      icel1_loc[nbr_fac_ini + ifac] = -1;

      /* icel2 */

      icel2_loc[nbr_fac_ini + ifac] = (ecs_int_t) *(int_lec + 1);

      /* nombres de sommets de la face */

      nbr_n_loc[nbr_fac_ini + ifac] = (ecs_int_t) *(int_lec + 5);

    }

    icoul_loc[nbr_fac_ini + ifac] = couleur;

    /* connectivite */

    cpt_som_total
      =  taille_connect_loc + cpt_som + nbr_n_loc[nbr_fac_ini + ifac];

    if (max_som < cpt_som_total) {
      ECS_REALLOC((*faces)->connect, 2 * max_som, ecs_int_t);
      max_som = 2 * max_som;
    }

    if (max_som_fac < nbr_n_loc[nbr_fac_ini + ifac]) {
      ECS_REALLOC(tab_lec, max_som_fac * 2, int32_t);
      max_som_fac = max_som_fac * 2;
    }

    ecs_file_read(tab_lec, 4, nbr_n_loc[nbr_fac_ini + ifac], fic);

    for (isom = 0; isom < nbr_n_loc[nbr_fac_ini + ifac]; isom ++) {
      (*faces)->connect[taille_connect_loc + cpt_som + isom]
        = (ecs_int_t) tab_lec[isom];
    }


    /* Mise a jour du compteur */

    cpt_som += nbr_n_loc[nbr_fac_ini + ifac];

  }


  /* Mise a jour de la taille de la connectivité */

  (*faces)->taille_connect = taille_connect_loc + cpt_som;
  ECS_REALLOC((*faces)->connect, taille_connect_loc + cpt_som, ecs_int_t);

  ECS_FREE(int_lec);
  ECS_FREE(tab_lec);

}


/*----------------------------------------------------------------------------
 *  Lecture d'une section cell et concatenation avec les eventuels
 *  elements appartenant aux sections cell prealablement lues
 *----------------------------------------------------------------------------*/

static void
ecs_pre_comet__lit_cell(ecs_file_t             *fic,
                        ecs_loc_cels_comet_t  **cels,
                        ecs_int_t              *taille_cels)
{
  ecs_int_t nbr_cel_loc;
  ecs_int_t nbr_cel_ini;
  ecs_int_t nbr_cel_total;
  ecs_int_t *icoul_loc;
  ecs_int_t *id_loc;

  ecs_int_t icel;

  int couleur;

  int32_t *int_lec = NULL;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  ECS_MALLOC(int_lec, 3, int32_t);

  ecs_file_read(int_lec, 4, 3, fic);
  nbr_cel_loc = (ecs_int_t ) *int_lec;
  couleur = *(int_lec + 1);

  printf(_("\n    Reading a cells section\n"));
  printf(_("         %10d cells "), (int)nbr_cel_loc);
  printf(_("of color %4d\n\n"), couleur);


  /* allocation */

  if (*taille_cels == 0) {
    *taille_cels = 1;
    ECS_REALLOC(*cels, *taille_cels, ecs_loc_cels_comet_t);
    nbr_cel_ini = 0;
    (*cels)->nbr_cel = nbr_cel_loc;
    ECS_MALLOC((*cels)->icoul, nbr_cel_loc, ecs_int_t);
    ECS_MALLOC((*cels)->id, nbr_cel_loc, ecs_int_t);
  }
  else {
    nbr_cel_ini = (*cels)->nbr_cel;
    nbr_cel_total = nbr_cel_ini + nbr_cel_loc;
    (*cels)->nbr_cel = nbr_cel_total;
    ECS_REALLOC((*cels)->icoul, nbr_cel_total, ecs_int_t);
    ECS_REALLOC((*cels)->id, nbr_cel_total, ecs_int_t);
  }

  id_loc = (*cels)->id;
  icoul_loc = (*cels)->icoul;

  ECS_REALLOC(int_lec, nbr_cel_loc, int32_t);

  ecs_file_read(int_lec, 4, nbr_cel_loc, fic);

  for (icel = 0; icel < nbr_cel_loc; icel ++) {
    id_loc[nbr_cel_ini + icel] = (ecs_int_t) int_lec[icel];
    icoul_loc[nbr_cel_ini + icel] = couleur;
  }

  ECS_FREE(int_lec);
}

/*----------------------------------------------------------------------------
 *  Lecture d'une section shell, les infos ne sont pas stockees
 *----------------------------------------------------------------------------*/

static void
ecs_pre_comet__lit_shell(ecs_file_t  *fic)
{
  ecs_int_t icel;

  int couleur;
  ecs_int_t count;

  int32_t *int_lec = NULL;
  ecs_int_t taille_lec;
  ecs_int_t nbr_shell_loc;

  /*Xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  taille_lec = 3;

  ECS_MALLOC(int_lec, taille_lec, int32_t);

  ecs_file_read(int_lec, 4, 3, fic);
  nbr_shell_loc = (ecs_int_t ) *int_lec;
  couleur = *(int_lec + 1);

  printf(_("\n    Reading a Shell section\n"));
  printf(_("\n    Caution: reading this type of section "
           "has not been validated\n"));

  for (icel = 0; icel < nbr_shell_loc; icel ++) {
    ecs_file_read(int_lec, 4, 2, fic);
    count = (ecs_int_t ) *(int_lec + 1);
    if (taille_lec < count) {
      ECS_MALLOC(int_lec, count, int32_t);
      taille_lec = count;
    }
    ecs_file_read(int_lec, 4, count, fic);
  }

  ECS_FREE(int_lec);
}

/*----------------------------------------------------------------------------
 *  Lecture d'une section map
 *----------------------------------------------------------------------------*/

static void
ecs_pre_comet__lit_map(ecs_file_t  *fic)
{
  ecs_int_t icel;

  int couleur;
  ecs_int_t tag;

  int32_t *int_lec = NULL;
  ecs_int_t taille_lec;
  ecs_int_t nbr_map_loc;

  char *chaine;

 /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  taille_lec = 3;

  ECS_MALLOC(int_lec, taille_lec, int32_t);

  ecs_file_read(int_lec, 4, 3, fic);
  nbr_map_loc = (ecs_int_t ) *int_lec;
  couleur = *(int_lec + 1);

  printf(_("\n    Reading a Map section\n"));
  printf(_("\n    Caution: reading this type of section "
           "has not been validated\n"));

  for (icel = 0; icel < nbr_map_loc; icel ++) {

    ecs_file_read(int_lec, 4, 1, fic);

    /* Lecture d'un entier */

    tag = (ecs_int_t ) *(int_lec + 1);

    /* Lecture d'une chaine */

    chaine = ecs_pre_comet__lit_chaine_xdr(fic);
    ECS_FREE(chaine);
  }

  ECS_FREE(int_lec);

}

/*----------------------------------------------------------------------------
 *  Lecture d'une section interface, les infos ne sont pas stockées
 *----------------------------------------------------------------------------*/

static void
ecs_pre_comet__lit_interf(ecs_file_t  *fic,
                          int          num_version)
{
  ecs_int_t icel;

  int couleur;

  int32_t *int_lec = NULL;
  ecs_int_t taille_lec;
  ecs_int_t nbr_shell_loc;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  taille_lec = 6;

  ECS_MALLOC(int_lec, taille_lec, int32_t);

  ecs_file_read(int_lec, 4, 3, fic);
  nbr_shell_loc = (ecs_int_t ) *int_lec;
  couleur = *(int_lec + 1);

  printf(_("\n    Reading an Interface section\n"));
  printf(_("\n    Caution: reading this type of section "
           "has not been validated\n"));

  for (icel = 0; icel < nbr_shell_loc; icel ++) {
    if (num_version >= 2)
      ecs_file_read(int_lec, 4, 6, fic);
    else
      ecs_file_read(int_lec, 4, 3, fic);
  }

  ECS_FREE(int_lec);
}

/*----------------------------------------------------------------------------
 *  Lecture d'une section scalar, les infos ne sont pas stockées
 *----------------------------------------------------------------------------*/

static void
ecs_pre_comet__lit_scalar(ecs_file_t  *fic)
{
  ecs_int_t icel;
  int couleur;

  int32_t *int_lec = NULL;
  ecs_int_t taille_lec;
  ecs_int_t nbr_shell_loc;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  taille_lec = 3;

  ECS_MALLOC(int_lec, taille_lec, int32_t);

  ecs_file_read(int_lec, 4, 3, fic);
  nbr_shell_loc = (ecs_int_t ) *int_lec;
  couleur = *(int_lec + 1);

  printf(_("\n    Reading a Scalar section\n"));
  printf(_("\n    Caution: reading this type of section "
           "has not been validated\n"));

  for (icel = 0; icel < nbr_shell_loc; icel ++)
    ecs_file_read(int_lec, 4, 2, fic);

  ECS_FREE(int_lec);
}

/*----------------------------------------------------------------------------
 * Extraction de la connectivité "cellules -> faces" d'un maillage.
 *
 * On considère une numérotation commune des faces internes et des
 * faces de bord, dans laquelle les faces de bord sont définies en
 * premier. L'indice commun de la i-ème face de bord est donc égal à i,
 * et celui de la j-ième face interne à nbr_fbr + j.
 *----------------------------------------------------------------------------*/

static void
ecs_pre_comet__cel_fac(const ecs_loc_faces_comet_t         *faces,
                       const ecs_loc_cels_comet_t          *cels,
                       ecs_int_t                    **const p_pos_cel_fac,
                       ecs_int_t                    **const p_val_cel_fac)
{
  ecs_int_t    icel, icel1, icel2, ifac, nbr_cel_loc;

  ecs_int_t  *cpt_cel_fac = NULL;
  ecs_int_t  *pos_cel_fac = NULL;
  ecs_int_t  *val_cel_fac = NULL;

  /* Allocation et initialisation de l'indice des positions */

  nbr_cel_loc = cels->nbr_cel;

  ECS_MALLOC(pos_cel_fac, nbr_cel_loc + 1, ecs_int_t);

  for (icel = 0; icel < nbr_cel_loc + 1; icel++)
    pos_cel_fac[icel] = 0;

  /* Comptage du nombre de faces par cellule
   * (on affecte le compteur temporaire correspondant à icel à
   * pos_cel_fac[icel + 1] et non pas à pos_cel_fac[icel] pour
   * faciliter l'étape suivante) */

  /* Remarque : test si icel < maillage->nbr_cel sur faces internes
     pour ignorer les cellules fantômes parallèles et/ou périodiques */

  for (ifac = 0; ifac < faces->nbr_fac; ifac++) {
    icel1 = faces->icel1[ifac] - 1;
    icel2 = faces->icel2[ifac] - 1;
    if (icel1 >= 0)
      pos_cel_fac[icel1 + 1] += 1;
    if (icel2 >= 0)
      pos_cel_fac[icel2 + 1] += 1;
  }

  /* Construction de l'indice des positions */

  pos_cel_fac[0] = 1;
  for (icel = 0; icel < nbr_cel_loc; icel++)
    pos_cel_fac[icel + 1] = pos_cel_fac[icel] + pos_cel_fac[icel + 1];

  /* Construction du tableau des valeurs */

  ECS_MALLOC(val_cel_fac, pos_cel_fac[nbr_cel_loc] - 1, ecs_int_t);
  ECS_MALLOC(cpt_cel_fac, nbr_cel_loc, ecs_int_t);

  for (icel = 0; icel < nbr_cel_loc; icel++)
    cpt_cel_fac[icel] = 0;

  /* Orientation des faces */

  /* Inversion de l'orientation pour les faces internes */

  for (ifac = 0; ifac < faces->nbr_fac; ifac++) {
    icel1 = faces->icel1[ifac] - 1;
    icel2 = faces->icel2[ifac] - 1;
    if (icel1 >= 0 ) {
      val_cel_fac[pos_cel_fac[icel1] + cpt_cel_fac[icel1] - 1]
        = (ifac + 1);
      cpt_cel_fac[icel1] += 1;
    }
    if (icel2 >= 0) {
      if  (icel1 >= 0 ) {
        val_cel_fac[pos_cel_fac[icel2] + cpt_cel_fac[icel2] - 1]
          =  -(ifac + 1);
        cpt_cel_fac[icel2] += 1;
      }
      else {
        val_cel_fac[pos_cel_fac[icel2] + cpt_cel_fac[icel2] - 1]
          =  ifac + 1;
        cpt_cel_fac[icel2] += 1;
      }
    }
  }
  ECS_FREE(cpt_cel_fac);

  /* Valeurs de retour */

  *p_pos_cel_fac = pos_cel_fac;
  *p_val_cel_fac = val_cel_fac;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  {
    ecs_int_t ipos, ival;
    /* Impression des tableaux */
    printf("dbg : cs_loc_maillage_ret_cel_fac\n"
           "nombre de cellules extraites = %d\n", nbr_cel_extr);
    for (ipos = 0; ipos < nbr_cel_extr; ipos++) {
      printf("  cellule %d\n", ipos);
      printf("    pos_cel_fac[%d] = %d\n", ipos, pos_cel_fac[ipos]);
      for (ival = pos_cel_fac[ipos]     - 1;
           ival < pos_cel_fac[ipos + 1] - 1;
           ival++)
        printf("      val_cel_fac[%d] = %d\n", ival, val_cel_fac[ival]);
    }
    printf("  pos_cel_fac[%d] = %d\n", ipos, pos_cel_fac[ipos]);
  }
#endif
}

/*----------------------------------------------------------------------------
 * Construction de la strucutre maillage
 *----------------------------------------------------------------------------*/

static ecs_maillage_t *
ecs_pre_comet__prepa_mail(ecs_loc_noeuds_comet_t  *noeuds,
                          ecs_loc_faces_comet_t   *faces,
                          ecs_loc_cels_comet_t    *cels,
                          ecs_int_t               *pos_cel_fac,
                          ecs_int_t               *val_cel_fac)
{
  /* Déclarations des variables de stockage        */
  /* avant transfert dans la structure du maillage */
  /*-----------------------------------------------*/

  size_t       cpt_elt_ent        [ECS_N_ENTMAIL]; /* Nbr. elts par entité */
  ecs_int_t    cpt_coul_ent       [ECS_N_ENTMAIL]; /* Compteur de couleurs */
  size_t       cpt_val_som_ent    [ECS_N_ENTMAIL]; /* Taille connect.      */
  ecs_int_t   *val_coul_ent       [ECS_N_ENTMAIL]; /* Tableau des couleurs */
  ecs_size_t  *cpt_elt_coul_ent   [ECS_N_ENTMAIL];
  ecs_size_t  *elt_pos_som_ent    [ECS_N_ENTMAIL]; /* Positions sommets    */
  ecs_int_t   *elt_val_som_ent    [ECS_N_ENTMAIL]; /* Numéros des sommets  */
  ecs_int_t   *elt_val_color_ent  [ECS_N_ENTMAIL]; /* Couleurs éléments    */

  ecs_int_t ifac;
  ecs_int_t nbr_som_elt;
  ecs_int_t icoul;

  ecs_entmail_t  entmail_e;

  ecs_int_t      ient;
  ecs_int_t      pos_elt;
  ecs_int_t      isom;
  ecs_int_t      nbr_fac;
  ecs_int_t      icel;
  ecs_int_t      pos_fac;
  ecs_int_t      nbr_som_fac;
  ecs_int_t      num_som_deb_fac;

  /* Création d'un maillage initialement vide (valeur de retour) */

  ecs_maillage_t  *maillage = ecs_maillage__cree_nodal();

  /*Xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Sommets */

  ecs_maillage_pre__cree_som(maillage,
                             noeuds->nbr_noeuds,
                             noeuds->coord);

  /* Elements */

  for (ient = ECS_ENTMAIL_FAC; ient < ECS_N_ENTMAIL; ient++) {

    cpt_elt_ent        [ient] = 0;
    cpt_val_som_ent    [ient] = 0;

    elt_pos_som_ent    [ient] = NULL;
    elt_val_som_ent    [ient] = NULL;
    elt_val_color_ent  [ient] = NULL;

    cpt_coul_ent       [ient] = 0;
    val_coul_ent       [ient] = NULL;
    cpt_elt_coul_ent   [ient] = NULL;

  }

  /* Faces (polygones) */

  entmail_e = ECS_ENTMAIL_FAC;

  cpt_elt_ent[entmail_e] = faces->nbr_fac;
  cpt_val_som_ent[entmail_e]  = faces->taille_connect;
  elt_val_som_ent[entmail_e]  = faces->connect;

  if (cpt_elt_ent[entmail_e] > 0) {

    ECS_MALLOC(elt_pos_som_ent[entmail_e],
               cpt_elt_ent[entmail_e] + 1,
               ecs_size_t);

    elt_pos_som_ent[entmail_e][0] = 1;

    ECS_MALLOC(elt_val_color_ent[entmail_e]  ,
               cpt_elt_ent[entmail_e],
               ecs_int_t);
  }

  for (ifac = 0; ifac < faces->nbr_fac; ifac  ++) {

    /* Affectation couleurs */

    for (icoul = 0;
            icoul < cpt_coul_ent[entmail_e]
         && val_coul_ent[entmail_e][icoul] != faces->icoul[ifac];
         icoul++);

    if (icoul == cpt_coul_ent[entmail_e]) {

      /* La valeur de la couleur n'a pas encore été stockée */

      ECS_REALLOC(val_coul_ent[entmail_e],
                  cpt_coul_ent[entmail_e] + 1,
                  ecs_int_t);
      ECS_REALLOC(cpt_elt_coul_ent[entmail_e],
                  cpt_coul_ent[entmail_e] + 1,
                  ecs_size_t);

      cpt_elt_coul_ent[entmail_e][icoul] = 0;
      val_coul_ent[entmail_e][icoul] = faces->icoul[ifac];
      cpt_coul_ent[entmail_e]++;

    }

    /* Affectation couleurs */

    cpt_elt_coul_ent[entmail_e][icoul]++;
    elt_val_color_ent[entmail_e][ifac] = icoul + 1;

    nbr_som_elt = faces->nbr_n[ifac];

    /* Construction connectivité */

    pos_elt = elt_pos_som_ent[entmail_e][ifac];

    elt_pos_som_ent[entmail_e][ifac + 1] = pos_elt + nbr_som_elt;

  }

  /* Cellules (polyèdres) */

  entmail_e = ECS_ENTMAIL_CEL;

  cpt_elt_ent[entmail_e] = cels->nbr_cel;

  nbr_fac = pos_cel_fac[cels->nbr_cel] - 1;

  cpt_val_som_ent[entmail_e] = 0;

  for (ifac = 0; ifac < nbr_fac; ifac++)
    cpt_val_som_ent[entmail_e] += faces->nbr_n[abs(val_cel_fac[ifac])-1] + 1;

  if (cpt_elt_ent[entmail_e] > 0) {

    ECS_MALLOC(elt_pos_som_ent[entmail_e],
               cpt_elt_ent[entmail_e] + 1,
               ecs_size_t);

    elt_pos_som_ent[entmail_e][0] = 1;

    ECS_MALLOC(elt_val_som_ent[entmail_e],
               cpt_val_som_ent[entmail_e],
               ecs_int_t);

    ECS_MALLOC(elt_val_color_ent[entmail_e],
               cpt_elt_ent[entmail_e],
               ecs_int_t);
  }

  for (icel = 0; icel < cels->nbr_cel; icel++){

    /* Position dans la connectivité */

    cpt_val_som_ent[entmail_e] = 0;

    for (ifac = pos_cel_fac[icel] - 1;
         ifac < pos_cel_fac[icel + 1] - 1;
         ifac++)
      cpt_val_som_ent[entmail_e]
        += faces->nbr_n[abs(val_cel_fac[ifac]) - 1] + 1;

    elt_pos_som_ent[entmail_e][icel + 1]
      = elt_pos_som_ent[entmail_e][icel] + cpt_val_som_ent[entmail_e];

    /* Affectation couleurs */

    for (icoul = 0;
            icoul < cpt_coul_ent[entmail_e]
         && val_coul_ent[entmail_e][icoul] != cels->icoul[icel];
         icoul++);

    if (icoul == cpt_coul_ent[entmail_e]) {

      /* La valeur de la couleur n'a pas encore été stockée */

      ECS_REALLOC(val_coul_ent[entmail_e],
                  cpt_coul_ent[entmail_e] + 1,
                  ecs_int_t);
      ECS_REALLOC(cpt_elt_coul_ent[entmail_e],
                  cpt_coul_ent[entmail_e] + 1,
                  ecs_size_t);

      cpt_elt_coul_ent[entmail_e][icoul] = 0;
      val_coul_ent[entmail_e][icoul] = cels->icoul[icel];
      cpt_coul_ent[entmail_e]++;

    }

    /* Affectation couleurs */

    cpt_elt_coul_ent[entmail_e][icoul]++;
    elt_val_color_ent[entmail_e][icel] = icoul + 1;
  }

  /* Connectivité */

  /* Boucle sur toutes les faces de val_cel_fac */

  cpt_val_som_ent[entmail_e] = 0;

  for (ifac = 0; ifac < nbr_fac; ifac++) {

    pos_fac = elt_pos_som_ent[ECS_ENTMAIL_FAC][abs(val_cel_fac[ifac])-1];

    nbr_som_fac = faces->nbr_n[abs(val_cel_fac[ifac])-1];

    /* Orientation */

    if (val_cel_fac[ifac] < 0 ) {

      num_som_deb_fac = faces->connect[pos_fac + (nbr_som_fac - 1) -1];

      for (isom = 0; isom < nbr_som_fac; isom++)
        elt_val_som_ent[entmail_e][cpt_val_som_ent[entmail_e] + isom]
          = faces->connect[pos_fac -1 + ((nbr_som_fac - 1) - isom)];

    }

    else {

      num_som_deb_fac = faces->connect[pos_fac - 1];

      for (isom = 0; isom < nbr_som_fac; isom++)
        elt_val_som_ent[entmail_e][cpt_val_som_ent[entmail_e] + isom]
          = faces->connect[pos_fac -1 + isom];
    }

    elt_val_som_ent[entmail_e][cpt_val_som_ent[entmail_e] + isom]
      = num_som_deb_fac;

    cpt_val_som_ent[entmail_e] += nbr_som_fac + 1;
  }

  /* Transfert des valeurs lues dans les structures d'entité de maillage */

  ecs_maillage_pre__cree_elt(maillage,
                             cpt_elt_ent,
                             elt_pos_som_ent,
                             elt_val_som_ent,
                             NULL,
                             elt_val_color_ent,
                             cpt_coul_ent,
                             val_coul_ent,
                             cpt_elt_coul_ent);

  ECS_FREE(pos_cel_fac);
  ECS_FREE(val_cel_fac);

  /* Renvoi de la structure de maillage */

  return maillage;
}

/*============================================================================
 *                             Fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Lecture d'un fichier au format Comet
 *   et affectation des données dans la structure de maillage
 *----------------------------------------------------------------------------*/

ecs_maillage_t *
ecs_pre_comet__lit_maillage(const char  *nom_fic)
{
  ecs_file_t *fic;

  char *en_tete;

  ecs_int_t section_inconnue;

  ecs_int_t taille_noeuds;
  ecs_int_t taille_faces_internes;
  ecs_int_t taille_faces_bord;
  ecs_int_t taille_cels;

  ecs_loc_noeuds_comet_t *liste_noeuds;
  ecs_loc_faces_comet_t  *liste_faces_bord;
  ecs_loc_faces_comet_t  *liste_faces_internes;
  ecs_loc_faces_comet_t  *liste_faces;
  ecs_loc_cels_comet_t   *liste_cels;

  int32_t * type = NULL;
  ecs_int_t nbr_elts_lus;

  ecs_int_t num_version;

  bool        is_little_endian;
  unsigned int_endian;

  ecs_tab_int_t tab_connect;
  ecs_tab_int_t tab_label;

  ecs_int_t taille_connect_bord;
  ecs_int_t taille_connect_interne;
  ecs_int_t taille_connect;
  ecs_int_t nbr_fabord;
  ecs_int_t nbr_faint;
  ecs_int_t nbr_fac;
  ecs_int_t isom;
  ecs_int_t ifac;

  ecs_int_t *pos_cel_fac;
  ecs_int_t *val_cel_fac;

  /* Création d'un maillage initialement vide (valeur de retour) */

  ecs_maillage_t  *maillage = ecs_maillage__cree_nodal();

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


  /* Test si le systeme est little/big endian */

  int_endian = 0;
  is_little_endian = false;
  *((char *) (&int_endian)) = '\1';

  if (int_endian == 1)
    is_little_endian = true;

  printf(_("\n\n"
           "Reading mesh from file in pro-STAR/STAR4 format\n"
           "----------------------\n"));

  printf(_("  \"ngeom\" file: %s\n"),
         nom_fic);


  /* Ouverture du fichier ngeom */

  fic = ecs_file_open(nom_fic,
                      ECS_FILE_MODE_READ,
                      ECS_FILE_TYPE_BINARY);

  ecs_file_set_big_endian(fic);


  /* Initialisation */

  section_inconnue = 0;

  taille_noeuds = 0;
  liste_noeuds = NULL;

  taille_faces_bord = 0;
  liste_faces_bord  = NULL;

  taille_faces_internes = 0;
  liste_faces_internes  = NULL;

  taille_cels = 0;
  liste_cels  = NULL;

  /* En-tete */

  en_tete = ecs_pre_comet__lit_chaine_xdr(fic);
  printf(_("\n  Header: %s\n"),en_tete);
  ECS_FREE(en_tete);

  /* numero de version */

  ECS_MALLOC(type, 1, int32_t);

  nbr_elts_lus = ecs_file_read_try(type, 4, 1, fic);

  printf(_("\n  Version number: %10d\n"),*type);
  num_version = (ecs_int_t) *type;

  /* Boucle sur les sections */

  nbr_elts_lus = ecs_file_read_try(type, 4, 1, fic);

  while (nbr_elts_lus == 1 && section_inconnue == 0) {

    switch(*type) {

    case VERTEX:
      ecs_pre_comet__lit_vertex(fic,
                                &liste_noeuds,
                                &taille_noeuds,
                                is_little_endian);
      break;

    case VERTEX_DOUBLE:
      ecs_pre_comet__lit_vertex_dbl(fic,
                                    &liste_noeuds,
                                    &taille_noeuds,
                                    is_little_endian);
      break;

    case INTERNAL_FACE:
      ecs_pre_comet__lit_face(fic,
                              &liste_faces_internes,
                              &taille_faces_internes,
                              num_version);
      break;

    case BOUNDARY:
      ecs_pre_comet__lit_boundary(fic,
                                  &liste_faces_bord,
                                  &taille_faces_bord,
                                  num_version);
      break;

    case CELL:
      ecs_pre_comet__lit_cell(fic,
                              &liste_cels,
                              &taille_cels);
      break;

    case SHELL:
      ecs_pre_comet__lit_shell(fic);
      break;

    case INTERFACE:
      ecs_pre_comet__lit_interf(fic, num_version);
      break;

    case SCALAR:
      ecs_pre_comet__lit_scalar(fic);
      break;

    case VERTEX_MAP:
      ecs_pre_comet__lit_map(fic);
      break;

    case SHELL_MAP:
      ecs_pre_comet__lit_map(fic);
      break;

    case BOUNDARY_MAP:
      ecs_pre_comet__lit_map(fic);
      break;

    case INTERNAL_FACE_MAP:
      ecs_pre_comet__lit_map(fic);
      break;

    case CELL_MAP:
      ecs_pre_comet__lit_map(fic);
      break;

    case -1:
      printf(_("\n       Finished reading the usable part of the file\n\n"));
      section_inconnue = 1;
      break;

    default:
      printf(_("\n\n       unknown section type %2d\n"
               "          Premature end of reading\n\n"), *type);
      section_inconnue = 1;
      break;

    }

    /* Prochaine section */

    nbr_elts_lus = ecs_file_read_try(type, 4, 1, fic);

  }


  /* Libération du fichier */

  ecs_file_free(fic);
  ECS_FREE(type);


  /* le passage en indices se fait avant la concaténation des
     faces pour conserver icel1 = -1 pour les faces de bord */

  /* labels des sommets */

  /* faces de bord */

  tab_connect.nbr = liste_faces_bord->taille_connect;
  tab_connect.val = liste_faces_bord->connect;

  tab_label.nbr = liste_noeuds->nbr_noeuds;
  tab_label.val = liste_noeuds->id_noeud;

  ecs_tab_int__ref_en_indice(tab_connect, tab_label, false);
  liste_faces_bord->connect   = tab_connect.val;

  /* faces internes */

  tab_connect.nbr = liste_faces_internes->taille_connect;
  tab_connect.val = liste_faces_internes->connect;

  ecs_tab_int__ref_en_indice(tab_connect, tab_label, false);
  liste_faces_internes->connect   = tab_connect.val;

  /* labels des cellules */

  tab_label.nbr = liste_cels->nbr_cel;
  tab_label.val = liste_cels->id;
  tab_connect.nbr = liste_faces_bord->nbr_fac;

  /* Pour les faces de bord, on ne traite que icel2 (icel1 = -1)*/

  tab_connect.val = liste_faces_bord->icel2;
  ecs_tab_int__ref_en_indice(tab_connect, tab_label, false);
  liste_faces_bord->icel2   = tab_connect.val;

  /* Faces internes */

  tab_connect.nbr = liste_faces_internes->nbr_fac;

  tab_connect.val = liste_faces_internes->icel1;
  ecs_tab_int__ref_en_indice(tab_connect, tab_label, false);
  liste_faces_internes->icel1   = tab_connect.val;

  tab_connect.val = liste_faces_internes->icel2;
  ecs_tab_int__ref_en_indice(tab_connect, tab_label, false);
  liste_faces_internes->icel2   = tab_connect.val;

  /* décalage 0 ... n-1 a 1 ... n */

  for(isom = 0; isom < liste_faces_internes->taille_connect; isom++)
    liste_faces_internes->connect[isom]
      = liste_faces_internes->connect[isom] + 1;

  for(isom = 0; isom < liste_faces_bord->taille_connect; isom++)
    liste_faces_bord->connect[isom] = liste_faces_bord->connect[isom] + 1;

  for(ifac = 0; ifac < liste_faces_internes->nbr_fac; ifac++){
    if (liste_faces_internes->icel2[ifac] != -1 )
      liste_faces_internes->icel2[ifac]
        = liste_faces_internes->icel2[ifac] + 1;
    if (liste_faces_internes->icel1[ifac] != -1 )
      liste_faces_internes->icel1[ifac]
        = liste_faces_internes->icel1[ifac] + 1;
  }

  for(ifac = 0; ifac < liste_faces_bord->nbr_fac; ifac++)
    if (liste_faces_bord->icel2[ifac] != -1 )
      liste_faces_bord->icel2[ifac] = liste_faces_bord->icel2[ifac] + 1;

  /* Concaténation des faces de bord et
     faces internes : faces de bord en premier */

  taille_connect_bord = liste_faces_bord->taille_connect;
  taille_connect_interne = liste_faces_internes->taille_connect;

  taille_connect = taille_connect_bord + taille_connect_interne;

  nbr_fabord = liste_faces_bord->nbr_fac;
  nbr_faint  = liste_faces_internes->nbr_fac;

  nbr_fac = nbr_fabord + nbr_faint;

  ECS_REALLOC(liste_faces_bord->connect, taille_connect, ecs_int_t);
  ECS_REALLOC(liste_faces_bord->icoul, nbr_fac, ecs_int_t);
  ECS_REALLOC(liste_faces_bord->icel1, nbr_fac, ecs_int_t);
  ECS_REALLOC(liste_faces_bord->icel2, nbr_fac, ecs_int_t);
  ECS_REALLOC(liste_faces_bord->nbr_n, nbr_fac, ecs_int_t);

  liste_faces = liste_faces_bord;

  for (isom = 0; isom < taille_connect_interne; isom++) {

    liste_faces->connect[taille_connect_bord + isom ]
      = liste_faces_internes->connect[isom];

  }

  for (ifac = 0; ifac < nbr_faint; ifac ++){

    liste_faces->icoul[nbr_fabord + ifac] = liste_faces_internes->icoul[ifac];
    liste_faces->icel1[nbr_fabord + ifac] = liste_faces_internes->icel1[ifac];
    liste_faces->icel2[nbr_fabord + ifac] = liste_faces_internes->icel2[ifac];
    liste_faces->nbr_n[nbr_fabord + ifac] = liste_faces_internes->nbr_n[ifac];

  }

  liste_faces->taille_connect = taille_connect;
  liste_faces->nbr_fac = nbr_fac;

  /* Libération de la memoire */

  ECS_FREE(liste_faces_internes->icoul);
  ECS_FREE(liste_faces_internes->icel1);
  ECS_FREE(liste_faces_internes->icel2);
  ECS_FREE(liste_faces_internes->nbr_n);
  ECS_FREE(liste_faces_internes->connect);

  ECS_FREE(liste_faces_internes);


  /* Préparation de la connectivité cellule */

  ecs_pre_comet__cel_fac(liste_faces,
                         liste_cels,
                         &pos_cel_fac,
                         &val_cel_fac);

  /* Remplissage des tableaux */

  maillage = ecs_pre_comet__prepa_mail(liste_noeuds,
                                       liste_faces,
                                       liste_cels,
                                       pos_cel_fac,
                                       val_cel_fac);

  ECS_FREE(liste_noeuds->id_noeud);
  ECS_FREE(liste_noeuds);

  ECS_FREE(liste_faces->icel1);
  ECS_FREE(liste_faces->icel2);
  ECS_FREE(liste_faces->nbr_n);
  ECS_FREE(liste_faces->icoul);

  ECS_FREE(liste_faces);

  ECS_FREE(liste_cels->icoul);
  ECS_FREE(liste_cels->id);

  ECS_FREE(liste_cels);

  /* Renvoi de la structure de maillage */

  return maillage;
}

/*----------------------------------------------------------------------------*/

