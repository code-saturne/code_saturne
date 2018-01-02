/*============================================================================
 *  Définition de la fonction
 *   de lecture d'un fichier de maillage NOPO (INRIA, utilisé par Simail)
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

/*============================================================================
 *                                 Visibilité
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie standard C
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdlib.h>
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
#include "ecs_maillage.h"
#include "ecs_maillage_priv.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_maillage_pre.h"


/*----------------------------------------------------------------------------
 *  Fichier  `include' du  paquetage courant associé au fichier courant
 *----------------------------------------------------------------------------*/

#include "ecs_pre_nopo.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' privés   du  paquetage courant
 *----------------------------------------------------------------------------*/


/*============================================================================
 *                  Définitions de paramètres et macros
 *============================================================================*/

#define ECS_NOPO_NBR_MAX_SOM           8
#define ECS_NOPO_NBR_MAX_SSELT         8
#define ECS_NOPO_NBR_MAX_TAB_REF      26        /* Valeur max de NMAE + 1
                                                   (nb. sommets + nb. arêtes
                                                   +nb. faces : 8 + 12 + 6) */

#define ECS_NOPO_NUL                   0        /* Inexistant */
#define ECS_NOPO_NODE                  1        /* Noeud */
#define ECS_NOPO_SEGMENT               2        /* Segment */
#define ECS_NOPO_TRIANGLE              3        /* Triangle */
#define ECS_NOPO_QUADRANGLE            4        /* Quadrangle */
#define ECS_NOPO_TETRAHEDRON           5        /* Tétraèdre */
#define ECS_NOPO_PENTAHEDRON           6        /* Prisme */
#define ECS_NOPO_HEXAHEDRON            7        /* Hexaèdre */
#define ECS_NOPO_SUPER_ELEMENT         8        /* Super-élément (inutilisé) */


/*============================================================================
 *                  Définition de structures locales
 *============================================================================*/

/* Définition des éléments */
/*=========================*/

typedef struct {

  ecs_int_t       nopo_typ;                       /* Type NOPO de l'élément  */
  ecs_elt_typ_t   ecs_typ;                        /* Type ECS  de l'élément  */
  ecs_int_t       num_som[ECS_NOPO_NBR_MAX_SOM];  /* Numéros de sommets ECS  */
  ecs_int_t       nbr_are_vol;                    /* Nombre arêtes si volume */
  ecs_int_t       nbr_sselt;                      /* Nombre de sous-éléments */
  ecs_sous_elt_t  sous_elt[ECS_NOPO_NBR_MAX_SSELT];

} ecs_loc_nopo_elt_t;

static const ecs_loc_nopo_elt_t  ecs_loc_nopo_elt_liste_c[8] = {

  {                        /* 1 */
    ECS_NOPO_NUL,
    ECS_ELT_TYP_NUL,
    { 0 },
    0,
    0,
    {
      {0,{0}}
    }
  },
  {                        /* 1 */
    ECS_NOPO_NODE,
    ECS_ELT_TYP_NUL,
    { 0 },
    0,
    0,
    {
      {0,{0}}
    }
  },
  {                        /* 2 */
    ECS_NOPO_SEGMENT,
    ECS_ELT_TYP_NUL,
    { 0 },
    0,
    0,
    {
      {0, {0}}                                     /*    1       2            */
    }                                              /*    x-------x            */
  },
  {                        /* 3 */
    ECS_NOPO_TRIANGLE,
    ECS_ELT_TYP_FAC_TRIA,
    { 1, 2, 3 },
    0,                                             /*        x 3              */
    0,                                             /*       / \               */
    {                                              /*      /   \              */
      {0, {0}}                                     /*     /     \             */
    }                                              /*  1 x-------x 2          */
  },
  {                        /* 4 */
    ECS_NOPO_QUADRANGLE,
    ECS_ELT_TYP_FAC_QUAD,
    { 1, 2, 3, 4 },
    0,                                             /*  4 x-------x 3          */
    0,                                             /*    |       |            */
    {                                              /*    |       |            */
      {0, {0}}                                     /*    |       |            */
    }                                              /*  1 x-------x 2          */
  },
  {                        /* 5 */
    ECS_NOPO_TETRAHEDRON,
    ECS_ELT_TYP_CEL_TETRA,
    { 1, 2, 3, 4 },
    6,                                             /*        x 4              */
    4,                                             /*       /|\               */
    {                                              /*      / | \              */
      {ECS_ELT_TYP_FAC_TRIA, { 1 , 3 , 2 }},       /*     /  |  \             */
      {ECS_ELT_TYP_FAC_TRIA, { 1 , 4 , 3 }},       /*  1 x- -|- -x 3          */
      {ECS_ELT_TYP_FAC_TRIA, { 1 , 2 , 4 }},       /*     \  |  /             */
      {ECS_ELT_TYP_FAC_TRIA, { 2 , 3 , 4 }}        /*      \ | /              */
    }                                              /*       \|/               */
  },                                               /*        x 2              */
  {                       /*  6 */
    ECS_NOPO_PENTAHEDRON,
    ECS_ELT_TYP_CEL_PRISM,
    { 1 , 2 , 3 , 4 , 5 , 6 },
    9,                                             /*  4 x-------x 6          */
    5,                                             /*    |\     /|            */
    {                                              /*    | \   / |            */
      {ECS_ELT_TYP_FAC_TRIA, { 1 , 3 , 2 }    },   /*  1 x- \-/ -x 3          */
      {ECS_ELT_TYP_FAC_QUAD, { 1 , 4 , 6 , 3 }},   /*     \ 5x  /             */
      {ECS_ELT_TYP_FAC_QUAD, { 1 , 2 , 5 , 4 }},   /*      \ | /              */
      {ECS_ELT_TYP_FAC_TRIA, { 4 , 5 , 6 }    },   /*       \|/               */
      {ECS_ELT_TYP_FAC_QUAD, { 2 , 3 , 6 , 5 }}    /*        x 2              */
    }
  },
  {                       /*  7 */
    ECS_NOPO_HEXAHEDRON,
    ECS_ELT_TYP_CEL_HEXA,
    { 1, 2, 3, 4, 5, 6, 7, 8 },
    12,
    6,                                             /*     8 x-------x 7       */
    {                                              /*      /|      /|         */
      {ECS_ELT_TYP_FAC_QUAD, { 1 , 4 , 3 , 2 }},   /*     / |     / |         */
      {ECS_ELT_TYP_FAC_QUAD, { 1 , 5 , 8 , 4 }},   /*  5 x-------x6 |         */
      {ECS_ELT_TYP_FAC_QUAD, { 1 , 2 , 6 , 5 }},   /*    | 4x----|--x 3       */
      {ECS_ELT_TYP_FAC_QUAD, { 5 , 6 , 7 , 8 }},   /*    | /     | /          */
      {ECS_ELT_TYP_FAC_QUAD, { 2 , 3 , 7 , 6 }},   /*    |/      |/           */
      {ECS_ELT_TYP_FAC_QUAD, { 3 , 4 , 8 , 7 }}    /*  1 x-------x 2          */
    }
  }
};

/*============================================================================
 *                              Fonctions privées
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Compactage des références. On fournit en entrée le numéro de référence
 *  d'une liste d'entités, on obtient en sortie (dans le même tableau)
 *  la position de cette référence dans le tableau compacté ; on récupère
 *  le tableau des valeurs des références correspondant à chaque indice,
 *  ainsi que le nombre de références total et le nombre d'entités par
 *  référence.
 *
 *  Si on n'a aucune référence, on renvoie des pointeurs à NULL.
 *----------------------------------------------------------------------------*/

static void
ecs_loc_pre_nopo__compct_ref(ecs_int_t     nbr_ent,
                             ecs_int_t    *cpt_coul_ent,
                             ecs_int_t   **elt_val_coul_ent,
                             ecs_size_t  **cpt_elt_coul_ent,
                             ecs_int_t   **val_coul_ent)
{
  ecs_int_t   ind;
  ecs_int_t   ind_cpt;
  ecs_int_t  *renum;

  int num_coul_min = 0;
  int num_coul_max = 0;
  int bw = 0;
  int ind_0 = -1;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (nbr_ent == 0)
    return;

  /* Détermination références min et max */

  num_coul_min = (*elt_val_coul_ent)[0];
  num_coul_max = (*elt_val_coul_ent)[0];

  for (ind = 1; ind < nbr_ent; ind++) {
    int elt_coul = (*elt_val_coul_ent)[ind];
    if (elt_coul > num_coul_max)
      num_coul_max = elt_coul;
    else if (elt_coul < num_coul_min)
      num_coul_min = elt_coul;
  }

  bw = num_coul_max - num_coul_min + 1;

  /* Indice correspondant à la référence 0 */

  if (num_coul_min <= 0)
    ind_0 = - num_coul_min;

  /* Compteurs */

  ECS_MALLOC(*cpt_elt_coul_ent, bw, ecs_size_t);

  for (ind = 0; ind < bw; ind++)
    (*cpt_elt_coul_ent)[ind] = 0;

  for (ind = 0; ind < nbr_ent; ind++)
    (*cpt_elt_coul_ent)[(*elt_val_coul_ent)[ind] - num_coul_min] += 1;

  /* Compactage des références réellement utilisées */

  ECS_MALLOC(*val_coul_ent, bw, ecs_int_t);
  ECS_MALLOC(renum, bw, ecs_int_t);

  ind_cpt = 0;

  for (ind = 0; ind < bw; ind++) {

    if ((*cpt_elt_coul_ent)[ind] != 0) {

      (*val_coul_ent)[ind_cpt] = ind + num_coul_min;
      (*cpt_elt_coul_ent)[ind_cpt] = (*cpt_elt_coul_ent)[ind];

      /* On ajoute 1 pour une numérotation de 1 à n et non de 0 à n-1 */

      renum[ind] = ind_cpt + 1;

      ind_cpt += 1;

    }
  }

  *cpt_coul_ent = ind_cpt;

  for (ind = 0; ind < nbr_ent; ind++)
    (*elt_val_coul_ent)[ind] = renum[(*elt_val_coul_ent)[ind] - num_coul_min];

  ECS_FREE(renum);

  ECS_REALLOC(*cpt_elt_coul_ent, *cpt_coul_ent, ecs_size_t);

  /* Test si on n'a que la référence 0 (équivaut à aucune référence) */

  if (*cpt_coul_ent == 1) {

    if ((*val_coul_ent)[0] == ind_0) {

      *cpt_coul_ent = 0;

      ECS_FREE(*cpt_elt_coul_ent);
      ECS_FREE(*elt_val_coul_ent);
      ECS_FREE(*val_coul_ent);

    }
  }
}

/*----------------------------------------------------------------------------
 *  Lecture de la table de connectivité ; on lit les éléments présents dans
 *  le maillage, mais on y ajoute des éléments correspondant aux
 *  sous-éléments référencés, afin de porter ces références. On ne traite
 *  ainsi que les sous-éléments de niveau directement inféreieur ;
 *  Ainsi, si deux faces et trois arêtes d'un tétraèdre donné sont référencées
 *  (colorées), on créera en plus du tétraèdre deux triangles portant
 *  les  références des faces, mais aucune arête.
 *----------------------------------------------------------------------------*/

static void
ecs_loc_pre_nopo__cree_elements(ecs_maillage_t   *maillage,
                                ecs_int_t        *nop5,
                                ecs_coord_t     **som_val_coord,
                                ecs_int_t         ncopnp,
                                ecs_int_t         ne,
                                ecs_int_t         np)
{
  ecs_int_t    def_som_elt[ECS_NOPO_NBR_MAX_SOM];

  ecs_int_t    nbr_som;
  ecs_int_t    nbr_sselt;

  ecs_int_t    nbr_pos_ent[ECS_N_ENTMAIL];
  ecs_int_t    nbr_val_ent[ECS_N_ENTMAIL];
  ecs_int_t    nbr_pos_add[ECS_N_ENTMAIL];
  ecs_int_t    nbr_val_add[ECS_N_ENTMAIL];

  ecs_int_t    cpt_coul_ent[ECS_N_ENTMAIL];
  ecs_int_t   *val_coul_ent[ECS_N_ENTMAIL];
  ecs_size_t  *cpt_elt_coul_ent[ECS_N_ENTMAIL];

  ecs_int_t    ient;
  ecs_int_t    ind;
  ecs_int_t    iel;
  ecs_int_t    iloc;
  ecs_int_t    ipos;
  ecs_int_t    isom;
  ecs_int_t    isselt;
  ecs_int_t    issent;

  ecs_int_t    tsselt;

  ecs_int_t    ncge;          /* Type d'élément */
  ecs_int_t    nmae;
  ecs_int_t    ndsde;
  ecs_int_t    nno;
  ecs_int_t    npo;
  ecs_int_t    ining;
  ecs_int_t    iref;
  ecs_int_t    tab_ref[ECS_NOPO_NBR_MAX_TAB_REF];

  /* Stockage avant transfert */
  /*--------------------------*/

  size_t       cpt_elt_ent        [ECS_N_ENTMAIL]; /* Nombre d'elems/entite */
  ecs_int_t    ind_elt_add        [ECS_N_ENTMAIL]; /* Indice ajouts/entité  */

  ecs_size_t  *elt_pos_som_ent    [ECS_N_ENTMAIL]; /* Positions numeros som */
  ecs_int_t   *elt_val_som_ent    [ECS_N_ENTMAIL]; /* Numeros des sommets   */
  ecs_int_t   *elt_val_coul_ent   [ECS_N_ENTMAIL]; /* Références            */

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /*====================================================*/
  /* Initialisations et allocations des tableaux locaux */
  /*====================================================*/

  for (ient = 0; ient < ECS_N_ENTMAIL; ient++) {

    cpt_elt_ent[ient] = 0;
    nbr_pos_ent[ient] = 1; /* Positions et valeurs éléments principaux */
    nbr_val_ent[ient] = 0;
    nbr_pos_add[ient] = 0; /* Positions et valeurs sous-éléments */
    nbr_val_add[ient] = 0;
  }

#define ECS_FCT_TYP(ncge) ecs_loc_nopo_elt_liste_c[ncge].ecs_typ

  /*========================================================*/
  /* Première boucle sur les éléments :                     */
  /* - comptages pour le dimensionnement des éléments       */
  /* - traitement des références des sommets                */
  /*========================================================*/

  ind = 0;

  for (iel = 0; iel < ne; iel++) {

    ncge   = nop5[ind++];

    /* Calcul de l'entité correspondante */

    if (ncge == ECS_NOPO_NODE)
      ient = ECS_ENTMAIL_NONE;
    else if (ncge == ECS_NOPO_SEGMENT)
      ient = ECS_N_ENTMAIL;
    else if (ncge < ECS_NOPO_TETRAHEDRON)
      ient = ECS_ENTMAIL_FAC;
    else if (ncge < ECS_NOPO_SUPER_ELEMENT)
      ient = ECS_ENTMAIL_CEL;
    else
      ient = ECS_N_ENTMAIL;

    nmae   = nop5[ind++]; /* Nombre de mots pour réf. faces, arêtes, ... */
    ndsde  = nop5[ind++]; /* Numéro de sous-domaine */
    nno    = nop5[ind++]; /* Nombre de noeuds pour l'élément */

    nbr_som = ecs_fic_elt_typ_liste_c[ECS_FCT_TYP(ncge)].nbr_som;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
    printf("ele %d : t %d; nmae %d; ndsde %d; nno %d (nbse %d)\n",
           iel + 1, ncge, nmae, ndsde, nno, nbr_som);
#endif

    for (isom = 0; isom < nno && isom < nbr_som; isom++)
      def_som_elt[isom] = nop5[ind++];
    for (        ; isom < nno                  ; isom++)
      ind++;

    if (ncopnp == 0) {     /* Si sommets non confondus avec les noeuds,
                              on prend les sommets (->éléments linéaires) */
      npo = nop5[ind++];
      assert (npo == nbr_som);
      for (isom = 0; isom < npo; isom++)
        def_som_elt[isom] = nop5[ind++];

    }

    if (nmae != 0) {

      ining = nop5[ind++];

      /* D'après la documention : boucle de 2 à NMAE (équiv. 0 à NMAE - 2) */

      for (iref = 0; iref < nmae - 1; iref++)
        tab_ref[iref] = nop5[ind++];

      /* En fonction de INING, on compte le nombre sous-éléments référencés */

      if (   (ining  < 3 && ient == ECS_ENTMAIL_FAC)
          || (ining == 1 && ient == ECS_ENTMAIL_CEL)) {

        nbr_sselt = ecs_loc_nopo_elt_liste_c[ncge].nbr_sselt;

        for (isselt = 0; isselt < nbr_sselt; isselt++) {

          if (tab_ref[isselt] != 0) {

            tsselt = ecs_loc_nopo_elt_liste_c[ncge].sous_elt[isselt].elt_typ;

            nbr_pos_add[ient - 1] += 1;
            nbr_val_add[ient - 1] += ecs_fic_elt_typ_liste_c[tsselt].nbr_som;

          }
        }
      }
    }

    /* Traitement selon l'entité correspondante */

    if (ient >= ECS_ENTMAIL_FAC && ient < ECS_N_ENTMAIL) {
      nbr_pos_ent[ient] += 1;
      nbr_val_ent[ient] += nbr_som;
    }

  }

  /*========================================================*/
  /* Création de la structure associée aux sommets          */
  /*========================================================*/

  /* Sommets */

  ecs_maillage_pre__cree_som(maillage, np, (*som_val_coord));

  /* Des tableaux sont libérés par ecs_maillage_pre__cree_som */

  *som_val_coord = NULL;

  /*========================================================*/
  /* Allocations restantes                                  */
  /*========================================================*/

  for (ient = ECS_ENTMAIL_FAC; ient < ECS_N_ENTMAIL; ient++) {

    if (nbr_pos_ent[ient] + nbr_pos_add[ient] > 1) {

      ECS_MALLOC(elt_pos_som_ent[ient],
                 nbr_pos_ent[ient] + nbr_pos_add[ient],
                 ecs_size_t);
      ECS_MALLOC(elt_val_som_ent[ient],
                 nbr_val_ent[ient] + nbr_val_add[ient],
                 ecs_int_t);
      ECS_MALLOC(elt_val_coul_ent[ient],
                 nbr_pos_ent[ient] + nbr_pos_add[ient] - 1,
                 ecs_int_t);

      elt_pos_som_ent[ient][0] = 1;

      /* Valeurs pos correspondant à fin éléments et début ajout */

      ind_elt_add[ient] = nbr_pos_ent[ient] - 1;

      elt_pos_som_ent[ient][ind_elt_add[ient]] = nbr_val_ent[ient] + 1;

    }
    else {

      elt_pos_som_ent[ient]     = NULL;
      elt_val_som_ent[ient]     = NULL;
      elt_val_coul_ent[ient]    = NULL;

      ind_elt_add[ient] = 0;

    }
  }

  /*========================================================*/
  /* Seconde boucle sur les éléments                        */
  /*========================================================*/

  ind = 0;

  for (iel = 0; iel < ne; iel++) {

    ncge   = nop5[ind++];

    /* Calcul de l'entité correspondante */

    if (ncge == ECS_NOPO_NODE)
      ient = ECS_ENTMAIL_NONE;
    else if (ncge == ECS_NOPO_SEGMENT)
      ient = ECS_N_ENTMAIL;
    else if (ncge < ECS_NOPO_TETRAHEDRON)
      ient = ECS_ENTMAIL_FAC;
    else if (ncge < ECS_NOPO_SUPER_ELEMENT)
      ient = ECS_ENTMAIL_CEL;
    else
      ient = ECS_N_ENTMAIL;

    nmae   = nop5[ind++]; /* Nombre de mots pour réf. faces, arêtes, ... */
    ndsde  = nop5[ind++]; /* Numéro de sous-domaine */
    nno    = nop5[ind++]; /* Nombre de noeuds pour l'élément */

    nbr_som = ecs_fic_elt_typ_liste_c[ECS_FCT_TYP(ncge)].nbr_som;

    for (isom = 0; isom < nno && isom < nbr_som; isom++)
      def_som_elt[isom] = nop5[ind++];
    for (        ; isom < nno                  ; isom++)
      ind++;

    if (ncopnp == 0) {     /* Si sommets non confondus avec les noeuds,
                              on prend les sommets (->éléments linéaires) */
      npo = nop5[ind++];
      for (isom = 0; isom < npo; isom++)
        def_som_elt[isom] = nop5[ind++];

    }

    if (nmae != 0) {

      ining = nop5[ind++];

      /* D'après la documention : boucle de 2 à NMAE (équiv. 0 à NMAE - 2) */

      for (iref = 0; iref < nmae - 1; iref++)
        tab_ref[iref] = nop5[ind++];

      /* En fonction de INING, on crée des sous-éléments référencés */

      if (   (ining  < 3 && ient == ECS_ENTMAIL_FAC)
          || (ining == 1 && ient == ECS_ENTMAIL_CEL)) {

        nbr_sselt = ecs_loc_nopo_elt_liste_c[ncge].nbr_sselt;

        for (isselt = 0; isselt < nbr_sselt; isselt++) {

          if (tab_ref[isselt] != 0) {

            /* Définition des sommets */

            tsselt = ecs_loc_nopo_elt_liste_c[ncge].sous_elt[isselt].elt_typ;

            issent = ient - 1;

            nbr_som = ecs_fic_elt_typ_liste_c[tsselt].nbr_som;

            ipos = elt_pos_som_ent[issent][ind_elt_add[issent]] - 1;

            for (isom = 0; isom < nbr_som; isom++) {

              iloc
                = ecs_loc_nopo_elt_liste_c[ncge].sous_elt[isselt].som[isom] - 1;

              elt_val_som_ent[issent][ipos + isom] = def_som_elt[iloc];

            }

            /* Position des numéros de sommets du prochain sous-élément */

            elt_pos_som_ent[issent][ind_elt_add[issent] + 1] =
              elt_pos_som_ent[issent][ind_elt_add[issent]] + nbr_som;

            /* Référence de l'élément ajouté (sera compactée plus tard) */

            elt_val_coul_ent[issent][ind_elt_add[issent]] = tab_ref[isselt];

            /* Incrémentation de l'indice d'éléments ajoutés */

            ind_elt_add[issent]++;

          }
        }
      }
    }

    /* Ajout de l'élément à l'entité correspondante (sauf pour les sommets) */

    if (ient >= ECS_ENTMAIL_FAC && ient < ECS_N_ENTMAIL) {

      nbr_som = ecs_fic_elt_typ_liste_c[ECS_FCT_TYP(ncge)].nbr_som;

      /* Connectivite de l'élément par ses numéros de sommets */

      for (isom = 0; isom <  nbr_som; isom++) {

        elt_val_som_ent
          [ient][elt_pos_som_ent[ient][cpt_elt_ent[ient]] - 1 + isom]
          = def_som_elt[ecs_loc_nopo_elt_liste_c[ncge].num_som[isom] - 1];
      }

      /* Position des numéros de sommets du prochain element */

      elt_pos_som_ent[ient][cpt_elt_ent[ient] + 1] =
        elt_pos_som_ent[ient][cpt_elt_ent[ient]] + nbr_som;

      /* Référence de l'élément lu (sera compactée plus tard) */

      elt_val_coul_ent[ient][cpt_elt_ent[ient]] = ndsde;

      /* Incrementation du nombre d'éléments lus */

      cpt_elt_ent[ient]++;

    }
  }

#undef ECS_FCT_TYP

  /* Mise à jour des dimensions */

  for (ient = ECS_ENTMAIL_FAC; ient < ECS_N_ENTMAIL; ient++)
    cpt_elt_ent[ient] = ind_elt_add[ient];

  /* Compactage des références des éléments */

  for (ient = ECS_ENTMAIL_FAC; ient < ECS_N_ENTMAIL; ient++) {

    if (cpt_elt_ent[ient] > 0) {

      ecs_loc_pre_nopo__compct_ref(cpt_elt_ent[ient],
                                   &(cpt_coul_ent[ient]),
                                   &(elt_val_coul_ent[ient]),
                                   &(cpt_elt_coul_ent[ient]),
                                   &(val_coul_ent[ient]));

    }
    else {

      ECS_FREE(elt_val_coul_ent[ient]);

      cpt_coul_ent[ient]     = 0;
      elt_val_coul_ent[ient] = NULL;
      cpt_elt_coul_ent[ient] = NULL;
      val_coul_ent[ient]     = NULL;
    }
  }

  /* Transfert des valeurs lues dans les structures d'entité de maillage */
  /*=====================================================================*/

  ecs_maillage_pre__cree_elt(maillage,
                             cpt_elt_ent,
                             elt_pos_som_ent,
                             elt_val_som_ent,
                             NULL,
                             elt_val_coul_ent,
                             cpt_coul_ent,
                             val_coul_ent,
                             cpt_elt_coul_ent);
}

/*----------------------------------------------------------------------------
 * Lecture d'un enregistrement d'un fichier NOPO. Chaque enregistrement
 * est un enregistrement de type binaire Fortran, constitué de mots de 4
 * ou 8 octets; dans le cas de mots de 4 octets, le premier mot donne la
 * dimension restante de l'enregistrement, et sera ignoré (i.e. pour un
 * tableau de n mots, on aura n + 1 entières, la première étant un entier
 * valant n).
 *
 * Le tableau des valeurs est alloué ici.
 *----------------------------------------------------------------------------*/

static void
ecs_loc_pre_nopo__lit_int(ecs_file_t     *fic_maillage,
                          const char     *nom_tableau,
                          size_t          elt_size,
                          size_t          n_elt,
                          size_t          n_seg,
                          size_t          l_seg,
                          ecs_int_t     **tableau)
{
  const char size_err_str[] = N_("Error reading NOPO file\n\"%s\"."
                                 "size of record \"%s\" is %d words\n"
                                 "but %d were expected.");

  ECS_MALLOC(*tableau, n_elt, ecs_int_t);

  if (elt_size == 4) {

    size_t i, j;
    size_t l_rem = n_elt - n_seg*l_seg;
    size_t l_max = ECS_MAX(l_seg, l_rem);
    int32_t *_tab = NULL;

    ECS_MALLOC(_tab, l_max + 1, int32_t);

    for (i = 0; i < n_seg; i++) {

      ecs_file_read(_tab, elt_size, l_seg + 1, fic_maillage);

      if ((size_t) (_tab[0]) != l_seg)
        ecs_error(__FILE__, __LINE__, 0, _(size_err_str),
                  ecs_file_get_name(fic_maillage), nom_tableau,
                  (int) ((*tableau)[0]), (int) l_seg);

      for (j = 0; j < l_seg; j++)
        (*tableau)[l_seg*i + j] = _tab[j+1];

    }

    if (l_rem > 0) {

      ecs_file_read(_tab, elt_size, l_rem + 1, fic_maillage);

      if ((size_t) (_tab[0]) != l_rem)
        ecs_error(__FILE__, __LINE__, 0, _(size_err_str),
                  ecs_file_get_name(fic_maillage), nom_tableau,
                  (int) ((*tableau)[0]), (int) l_rem);

      for (j = 0; j < l_rem; j++)
        (*tableau)[l_seg*n_seg + j] = _tab[j+1];

    }

    ECS_FREE(_tab);
  }

  else if (elt_size == 8) {

    size_t i, j;
    size_t l_rem = n_elt - n_seg*l_seg;
    size_t l_max = ECS_MAX(l_seg, l_rem);
    int64_t *_tab = NULL;

    ECS_MALLOC(_tab, l_max + 1, int64_t);

    for (i = 0; i < n_seg; i++) {

      ecs_file_read(_tab, elt_size, l_seg + 1, fic_maillage);

      if ((size_t) (_tab[0]) != l_seg)
        ecs_error(__FILE__, __LINE__, 0, _(size_err_str),
                  ecs_file_get_name(fic_maillage), nom_tableau,
                  (int) ((*tableau)[0]), (int) l_seg);

      for (j = 0; j < l_seg; j++)
        (*tableau)[l_seg*i + j] = _tab[j+1];
    }

    if (l_rem > 0) {

      ecs_file_read(_tab, elt_size, l_rem + 1, fic_maillage);

      if ((size_t) (_tab[0]) != l_rem)
        ecs_error(__FILE__, __LINE__, 0, _(size_err_str),
                  ecs_file_get_name(fic_maillage), nom_tableau,
                  (int) ((*tableau)[0]), (int) l_rem);

      for (j = 0; j < l_rem; j++)
        (*tableau)[l_seg*n_seg + j] = _tab[j+1];

    }

    ECS_FREE(_tab);
  }
}

static void
ecs_loc_pre_nopo__lit_real(ecs_file_t     *fic_maillage,
                           const char     *nom_tableau,
                           size_t          elt_size,
                           size_t          n_elt,
                           size_t          n_seg,
                           size_t          l_seg,
                           ecs_coord_t   **tableau)
{
  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  const char size_err_str[] = N_("Error reading NOPO file\n\"%s\"."
                                 "size of record \"%s\" is %d words\n"
                                 "but %d were expected.");

  ECS_MALLOC(*tableau, n_elt, ecs_coord_t);

  if (elt_size == 4) {

    size_t i, j;
    size_t l_rem = n_elt - n_seg*l_seg;
    size_t l_max = ECS_MAX(l_seg, l_rem);
    float *_tab = NULL;

    if (sizeof(float) != 4)
      ecs_error(__FILE__, __LINE__, 0,
                _("This function has been compiled with a "
                  "\"float\" type of size other than 4.\n"
                  "(Porting error making conversion of single-precision\n"
                  "NOPO coordinates to standard Code_Saturne Preprocessor\n"
                  "coordinates impossible)."));

    ECS_MALLOC(_tab, l_max + 1, float);

    for (i = 0; i < n_seg; i++) {

      ecs_file_read(_tab, elt_size, l_seg + 1, fic_maillage);

      if ((size_t) (((int32_t *)_tab)[0]) != l_seg)
        ecs_error(__FILE__, __LINE__, 0, _(size_err_str),
                  ecs_file_get_name(fic_maillage), nom_tableau,
                  (int) (((int32_t *)_tab)[0]), (int) l_seg);

      for (j = 0; j < l_seg; j++)
        (*tableau)[l_seg*i + j] = _tab[j+1];

    }

    if (l_rem > 0) {

      ecs_file_read(_tab, elt_size, l_rem + 1, fic_maillage);

      if ((size_t) (((int32_t *)_tab)[0]) != l_rem)
        ecs_error(__FILE__, __LINE__, 0, _(size_err_str),
                  ecs_file_get_name(fic_maillage), nom_tableau,
                  (int) (((int32_t *)_tab)[0]), (int) l_rem);

      for (j = 0; j < l_rem; j++)
        (*tableau)[l_seg*n_seg + j] = _tab[j+1];

    }

    ECS_FREE(_tab);
  }

  else if (elt_size == 8) {

    size_t i, j;
    size_t l_rem = n_elt - n_seg*l_seg;
    size_t l_max = ECS_MAX(l_seg, l_rem);
    double *_tab = NULL;

    if (sizeof(double) != 8)
      ecs_error(__FILE__, __LINE__, 0,
                _("This function has been compiled with a "
                  "\"double\" type of size other than 8.\n"
                  "(Porting error making conversion of double-precision\n"
                  "NOPO coordinates to standard Code_Saturne Preprocessor\n"
                  "coordinates impossible)."));

    ECS_MALLOC(_tab, l_max + 1, double);

    for (i = 0; i < n_seg; i++) {

      ecs_file_read(_tab, elt_size, l_seg + 1, fic_maillage);

      if ((size_t) (((int64_t *)_tab)[0]) != l_seg)
        ecs_error(__FILE__, __LINE__, 0, _(size_err_str),
                  ecs_file_get_name(fic_maillage), nom_tableau,
                  (int) (((int64_t *)_tab)[0]), (int) l_seg);

      for (j = 0; j < l_seg; j++)
        (*tableau)[l_seg*i + j] = _tab[j+1];

    }

    if (l_rem > 0) {

      ecs_file_read(_tab, elt_size, l_rem + 1, fic_maillage);

      if ((size_t) (((int64_t *)_tab)[0]) != l_rem)
        ecs_error(__FILE__, __LINE__, 0, _(size_err_str),
                  ecs_file_get_name(fic_maillage), nom_tableau,
                  (int) (((int64_t *)_tab)[0]), (int) l_rem);

      for (j = 0; j < l_rem; j++)
        (*tableau)[l_seg*n_seg + j] = _tab[j+1];

    }

    ECS_FREE(_tab);

  }
}

/*============================================================================
 *                             Fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Lecture d'un fichier NOPO (Format INRIA utilisé par Simail)
 *   et affectation des donnees dans la structure de maillage
 *----------------------------------------------------------------------------*/

ecs_maillage_t *
ecs_pre_nopo__lit_maillage(const char  *nom_fic_maillage)
{
  ecs_file_t    *fic_maillage;     /* Descripteur du fichier sur les noeuds */
  int            dim_e;            /* Dimension spatiale */

  ecs_coord_t  *som_val_coord;

  ecs_int_t     ind;
  int32_t       ind_test;

  size_t        lnop0 = 0; /* number of words in NOP0 */
  size_t        lnop1 = 0; /* number of words in NOP1 */
  size_t        lnop2 = 0; /* number of words in NOP2 */
  size_t        lnop3 = 0; /* number of words in NOP3 */
  size_t        lnop4 = 0; /* number of words in NOP4 */
  size_t        lnop5 = 0; /* number of words in NOP5 */

  ecs_int_t     *nop0;
  ecs_int_t     *nop1;
  ecs_int_t     *nop3;
  ecs_int_t     *nop2;
  ecs_coord_t   *nop4;
  ecs_int_t     *nop5;

  ecs_int_t     ndsr;    /* maximum reference number */
  ecs_int_t     ndsd;    /* maximum sub-domain number */
  ecs_int_t     ncopnp;  /* 1 if vertices = noeuds, 0 otherwise */
  ecs_int_t     ne;      /* number of elements */
  ecs_int_t     nepo;    /* number of point elements */
  ecs_int_t     nesg;    /* number of segments */
  ecs_int_t     ntri;    /* number of triangles */
  ecs_int_t     nqua;    /* number of quadrangles */
  ecs_int_t     ntet;    /* number of tetrahedra */
  ecs_int_t     npen;    /* number of pentahedra */
  ecs_int_t     nhex;    /* number of hexahedra */
  ecs_int_t     noe;     /* number of noeuds */
  ecs_int_t     np;      /* number of points */
  ecs_int_t     ntacoo;  /* coordinate system type:
                            1: cartesian; 2: cyl; 3: spherical */

  size_t  taille_e   = 56; /* 56 for 32-bit mode, 208 for 64-bit */
  size_t  taille_elt =  4; /* 4:32 bit or 6: 64 bit */

  size_t  ne0 = 0, le0 = 0; /* Segmentation of NOP0 array */
  size_t  ne1 = 0, le1 = 0; /* Segmentation of NOP1 array */
  size_t  ne2 = 0, le2 = 0; /* Segmentation of NOP2 array */
  size_t  ne3 = 0, le3 = 0; /* Segmentation of NOP3 array */
  size_t  ne4 = 0, le4 = 0; /* Segmentation of NOP4 array */
  size_t  ne5 = 0, le5 = 0; /* Segmentation of NOP5 array */

  int32_t nop0_c[29];

  /* Création d'un maillage initialement vide (valeur de retour) */

  ecs_maillage_t  *maillage = ecs_maillage__cree_nodal();

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Print title */
  /*=============*/

  printf(_("\n\n"
           "Reading mesh from file in NOPO (Simail) format\n"
           "----------------------\n"));

  printf(_("  Mesh file: %s\n\n\n"),
         nom_fic_maillage);

  /* Initialization */
  /*================*/

  /* Open NOPO file for reading */
  /*----------------------------*/

  fic_maillage = ecs_file_open(nom_fic_maillage,
                               ECS_FILE_MODE_READ,
                               ECS_FILE_TYPE_BINARY);

  /* Test if file is in native format or not */

  ecs_file_read((ecs_byte_t *)(&ind_test), sizeof(int32_t), 1,
                fic_maillage);

  if (ind_test != 56 && ind_test != 208) {

    ecs_file_swap_endian((ecs_byte_t *)(&ind_test), (ecs_byte_t *)(&ind_test),
                         sizeof(int32_t), 1);

    if (ind_test == 56 || ind_test == 208) {

      ecs_file_set_swap_endian
        (fic_maillage,
         (ecs_file_get_swap_endian(fic_maillage) == 1) ? 0 : 1);
    }
    else

      ecs_error(__FILE__, __LINE__, 0,
                _("Format error for file \"%s\" :\n"
                  "This file does not seem to be in NOPO format\n"
                  "(the header length is not 56 (32-bit) or 208 (64-bit))."),
                ecs_file_get_name(fic_maillage));
  }

  taille_e = ind_test;

  /* Read file header (array NOP0), defining the sizes of other arrays */

  /*
    Note that the first value of a NOP* array contains the array's size,
    so the array is thus of size n+1; We thus use 1 to n indexing
    rather than 0 to n-1 indexing to access other values.
  */

  if (taille_e == 56) { /* 32-bit */

    int32_t nopo_e[14];


    ecs_file_read((ecs_byte_t *)nopo_e, sizeof(int32_t),
                  taille_e / sizeof(int32_t),
                  fic_maillage);

    ecs_file_read((ecs_byte_t *)(&ind_test), sizeof(int32_t), 1,
                  fic_maillage);

    assert(ind_test == (int32_t)taille_e);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
    for (ind = 0; ind < 14; ind++)
      printf("nopo_e[%d] = %d\n", ind, nopo_e[ind]);
#endif

    /*
      According to the documentation, the header's 5th value is reserved
      and is 0, and the NOP3 array is unused. In practice, we note that
      the NOP3 array does exist when we have super-elements or descriptions,
      and the header's 5th value defines its size.
    */

    lnop0 = nopo_e[2];
    lnop1 = nopo_e[3];
    lnop2 = nopo_e[4];
    lnop3 = nopo_e[5];
    lnop4 = nopo_e[6];
    lnop5 = nopo_e[7];

    if (   nopo_e[ 1] !=  6 || nopo_e[ 2] != 32 || nopo_e[ 4] != 27
        || nopo_e[ 8] !=  1 || nopo_e[ 9] !=  1 || nopo_e[10] !=  1
        || nopo_e[11] !=  1 || nopo_e[12] !=  2 || nopo_e[13] !=  1)

      ecs_error(__FILE__, __LINE__, 0,
                _("Format error for file \"%s\" :\n"
                  "The reserved descriptor values are not those expected;\n"
                  "this file is probably not a NOPO/Simail file."),
                ecs_file_get_name(fic_maillage));

  }
  else if (taille_e == 208) { /* 64-bit */

    int64_t nopo_e[26];
    int32_t nopo_e_0[2];

    ecs_file_read(nopo_e, sizeof(int64_t),
                  taille_e / sizeof(int64_t),
                  fic_maillage);

    ecs_file_read((ecs_byte_t *)(&ind_test), sizeof(int32_t), 1,
                  fic_maillage);

    assert(ind_test == (int32_t)taille_e);

    /* The first 64-bit "integer" in fact contains 2 32-bit integers */

    if (ecs_file_get_swap_endian(fic_maillage) == 1) {
      ecs_file_swap_endian(&(nopo_e[0]), &(nopo_e[0]), 8, 1);
      ecs_file_swap_endian(&(nopo_e[0]), &(nopo_e[0]), 4, 2);
    }
    memcpy(nopo_e_0, &(nopo_e[0]), 8);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
    for (ind = 1; ind < 26; ind++)
      printf("nopo_e[%d] = %d\n", ind, nopo_e[ind]);
#endif

    /*
      The header's 5th value is usually 0. It corresponds to the
      size of array NOP3, which exists when we have super-elements
      or descriptions.
    */

    lnop0 = nopo_e[2];
    lnop1 = nopo_e[3];
    lnop2 = nopo_e[4];
    lnop3 = nopo_e[5];
    lnop4 = nopo_e[6];
    lnop5 = nopo_e[7];

    /* Array segmentation */

    ne0 = nopo_e[14]; le0 = nopo_e[15];
    ne1 = nopo_e[16]; le1 = nopo_e[17];
    ne2 = nopo_e[18]; le2 = nopo_e[19];
    ne3 = nopo_e[20]; le3 = nopo_e[21];
    ne4 = nopo_e[22]; le4 = nopo_e[23];
    ne5 = nopo_e[24]; le5 = nopo_e[25];

    /* According to the documentation, this array contains 26 64-bit
       records for 27 values: in fact, the 2 first values correspond
       to 32-bit values, encoded in a 64 bit integer. */

    if (nopo_e_0[1] == 64)
      taille_elt = 8;

    if (   nopo_e_0[0] != 26 || ((nopo_e_0[1] != 32) && (nopo_e_0[1] != 64))
        || nopo_e[ 1] !=  6 || nopo_e[ 2] != 32 || nopo_e[ 4] != 27
        || nopo_e[ 8] !=  1 || nopo_e[ 9] !=  1 || nopo_e[10] !=  1
        || nopo_e[11] !=  1 || nopo_e[12] !=  2 || nopo_e[13] !=  1)

      ecs_error(__FILE__, __LINE__, 0,
                _("Format error for file \"%s\" :\n"
                  "The reserved descriptor values are not those expected;\n"
                  "this file is probably not a NOPO/Simail file."),
                ecs_file_get_name(fic_maillage));

  }

  /* We may switch to Fortran binary once this record is read */

  ecs_file_set_type(fic_maillage, ECS_FILE_TYPE_FORTRAN_BINARY);

  /* Read NOP0 array */
  /* --------------- */

  ecs_loc_pre_nopo__lit_int(fic_maillage, "NOP0", taille_elt,
                            lnop0, ne0, le0, &nop0);

  /*
    NOP0 contains character strings encoded in integer arrays,
    so these should be restored if bytes were swapped.
    It seems that strings are encoded as 4-byte integers, so restoring them
    on a big-endian architecture for a file generated on a little-endian
    architecture requires the following call (the converse should be tested).
  */

  for (ind = 0; ind < 29; ind++)
    nop0_c[ind] = nop0[ind];

  if (ecs_file_get_swap_endian(fic_maillage) == 1) {
    ecs_file_swap_endian((ecs_byte_t *)nop0_c,
                         (ecs_byte_t *)nop0_c,
                         4, 29);
  }

  /* Remove whitespace at end of title (80 characters) */
  *((char *)(nop0_c) + 79) = '\0';
  for (ind = 79; ind > 0 && *((char *)(nop0_c) + ind) == ' '; ind--) {
    *((char *)(nop0_c) + ind) = '\0';
  }

  printf(_("  Title     : %.80s\n"), (char *)nop0_c);
  printf(_("  Date      : %2.2s/%2.2s/%4.4s\n"),
         (char *)(nop0_c + 20),
         ((char *)(nop0_c + 20)) + 2, ((char *)(nop0_c + 20)) + 4);
  printf(_("  Creator   : %24.24s\n"), (char *)(nop0_c + 22));

  if (strncmp("NOPO", (char *)(nop0_c + 28), 4) != 0)
    ecs_error(__FILE__, __LINE__, 0,
              _("Format error for file \"%s\" :\n"
                "String 'NOPO' does not appear in header;\n"
                "this file is probably not a NOPO/Simail file."),
              ecs_file_get_name(fic_maillage));

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  printf("nop0[   30] NIVEAU = %d\n",  nop0[29]);
  printf("nop0[   31] ETAT   = %d\n",  nop0[30]);
  printf("nop0[   32] NTACM  = %d\n",  nop0[31]);
#endif

  ECS_FREE(nop0);

  /* Read NOP1 array */
  /* --------------- */

  /*
   * This array contains axiliary arrays, not usually needed.
   */

  if (lnop1 != 0) {

    ecs_loc_pre_nopo__lit_int(fic_maillage, "NOP1",  taille_elt,
                              lnop1, ne1, le1, &nop1);

    ECS_FREE(nop1);
  }

  nop1    = NULL;

  /* Read NOP2 array */
  /* --------------- */

  ecs_loc_pre_nopo__lit_int(fic_maillage, "NOP2", taille_elt,
                            lnop2, ne2, le2, &nop2);

  /* Dimension du maillage */

  if (nop2[0] == 2)
    dim_e = 2;
  else
    dim_e = 3;

  printf(_("  Type      : %d bit\n"
           "  Dimension : %d\n\n"),
         (int)(taille_elt*8), (int) (nop2[0]));

  assert(nop2[0] == 2 || nop2[0] == 3);

  /* Other dimensions and parameters */

  ndsr   = (ecs_int_t) nop2[1];  /* maximum reference number */
  ndsd   = (ecs_int_t) nop2[2];  /* maximum sub-domain number */
  ncopnp = (ecs_int_t) nop2[3];  /* 1 if vertices = nodes, 0 otherwise */
  ne     = (ecs_int_t) nop2[4];  /* number of elements */
  nepo   = (ecs_int_t) nop2[5];  /* number of point elements */
  nesg   = (ecs_int_t) nop2[6];  /* nombre of segments */
  ntri   = (ecs_int_t) nop2[7];  /* nombre of triangles */
  nqua   = (ecs_int_t) nop2[8];  /* nombre of quadrangles */
  ntet   = (ecs_int_t) nop2[9];  /* nombre of tetrahedra */
  npen   = (ecs_int_t) nop2[10]; /* nombre of pentahedra */
  nhex   = (ecs_int_t) nop2[11]; /* number of hexahedra */
  noe    = (ecs_int_t) nop2[14]; /* number of noeuds */

  /*
   * Values not used here:
   *
   * Remark: according to the documentation, the number of words in NOP5
   *         is given by NOP2's 27th value, and by the header NOP0's 8th
   *         value. On some test cases (generated by Simail), the value
   *         given by NOP2 is not coherent; we ignore it, and keep the value
   *         read initially. In the same manner, the presence of NOP3
   *         depends on LNOP3 (5th value of NOP0), and not on NBEGM.
   *
   * nsup   = (ecs_int_t) nop2[12];  number of super-elements
   * nef    = (ecs_int_t) nop2[13];  number of boundary elements
   * n1     = (ecs_int_t) nop2[15];  n. interior nodes for segment or edge
   * iset   = (ecs_int_t) nop2[16];  n. interior nodes triangle ou face
   * iseq   = (ecs_int_t) nop2[17];  n. interior nodes for quadrangle or face
   * isete  = (ecs_int_t) nop2[18];  n. interior nodes for tetrahedra
   * isepe  = (ecs_int_t) nop2[19];  n. interior nodes for pentahedra
   * isehe  = (ecs_int_t) nop2[20];  n. interior nodes for hexahedra
   * ntycoo = (ecs_int_t) nop2[22];  coordinate value type (2 here)
   * lpgdn  = (ecs_int_t) nop2[23];  largest diff. between nodes of same elt.
   * nbegm  = (ecs_int_t) nop2[24];  number of super-elements in NOP3
   * lnop5  = (ecs_int_t) nop2[25];  number of words in NOP5
   */

  np     = (ecs_int_t) nop2[21];  /* number of points */

  ntacoo = (ecs_int_t) nop2[26];  /* coordinate system type:
                                     1: cartesian; 2: cyl; 3: spherical */

  ECS_FREE(nop2);

  printf(_("  Initial data: %10d points\n"
           "                %10d nodes\n"
           "                %10d point elements\n"
           "                %10d segments\n"
           "                %10d triangles\n"
           "                %10d quadrangles\n"),
         (int)np, (int)noe, (int)nepo, (int)nesg, (int)ntri, (int)nqua);

  if (dim_e == 2)
    printf("\n");
  else
    printf(_("                %10d tetrahedra\n"
             "                %10d pentahedra\n"
             "                %10d hexahedra\n\n"),
           (int)ntet, (int)npen, (int)nhex);

  if (ntacoo != 1)
    ecs_error(__FILE__, __LINE__, 0,
              _("Error reading NOPO file:\n\"%s\";\n"
                "The coordinate system is cylindrical or spherical,\n"
                "and this case is not currently handled"),
              ecs_file_get_name(fic_maillage));

  /* "Hardening" (maybe not necessary, but we prefer to be cautious) */

  assert(np <= noe);
  if (ncopnp == 1 && np == 0)
    np = noe;
  assert(np > 0);

  /* Read NOP3 array */
  /* --------------- */

  /* This array contains subarrays relative to super-elements */

  if (lnop3 != 0) {

    ecs_loc_pre_nopo__lit_int(fic_maillage, "NOP3",  taille_elt,
                              lnop3, ne3, le3, &nop3);

    ECS_FREE(nop3);
  }

  nop3 = NULL;

  /* Read NOP4 array */
  /* --------------- */

  ecs_loc_pre_nopo__lit_real(fic_maillage, "NOP4", taille_elt,
                             lnop4, ne4, le4, &nop4);

  /*
    "Z" coordinates may be absent from this array in 2D. as we are only
    interested by points and not by interior nodes, we only keep this
    part of the array (nodes come after).
  */

  if (dim_e == 3) {

    som_val_coord = nop4;
    nop4 = NULL;

  }
  else  if (dim_e == 2) {

    ECS_MALLOC(som_val_coord, np * 3, ecs_coord_t);

    for (ind = 0; ind < np; ind++) {
      som_val_coord[ind*3    ] = nop4[ind*2];
      som_val_coord[ind*3 + 1] = nop4[ind*2 + 1];
      som_val_coord[ind*3 + 2] = 0.0;
    }

    ECS_FREE(nop4);

  }

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  printf("Coordonnées\n");
  for (ind = 0; ind < np; ind++)
    printf("%d : % 10.5e % 10.5e % 10.5e\n",
           ind + 1, som_val_coord[ind*3    ],
           som_val_coord[ind*3 + 1], som_val_coord[ind*3 + 2]);
#endif

  /* Read NOP5 array */
  /* --------------- */

  ecs_loc_pre_nopo__lit_int(fic_maillage, "NOP5", taille_elt,
                            lnop5, ne5, le5, &nop5);

  /* Close mesh file */
  /*-----------------*/

  ecs_file_free(fic_maillage);

  /* Decode NOP5 array */
  /* ----------------- */

  ecs_loc_pre_nopo__cree_elements(maillage,
                                  nop5,
                                  &som_val_coord,
                                  ncopnp,
                                  ne,
                                  np);

  ECS_FREE(nop5);

  /* Return */
  /*--------*/

  return maillage;
}

/*----------------------------------------------------------------------------*/

