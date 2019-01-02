/*============================================================================
 *  Définition de la fonction
 *   de lecture d'un fichier de maillage IDEAS-MS au format "universel"
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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
 *  Fichiers `include' librairie standard C
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
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
#include "ecs_table.h"
#include "ecs_table_att.h"
#include "ecs_maillage.h"
#include "ecs_maillage_priv.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_maillage_pre.h"


/*----------------------------------------------------------------------------
 *  Fichier  `include' du  paquetage courant associé au fichier courant
 *----------------------------------------------------------------------------*/

#include "ecs_pre_ideas.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' privés   du  paquetage courant
 *----------------------------------------------------------------------------*/


/*============================================================================
 *  Définitions de parametres-macros
 *============================================================================*/

/* Pour une lecture de 80 caracteres par ligne  */
/* auxquels il faut ajouter le `\n' et le `\0'  */
/* pour l'affectation dans la chaine receptrice */
#define ECS_LOC_LNG_MAX_CHAINE_IDEAS  84        /* Dimension des chaines */

/* Numeros de dataset */
/*--------------------*/

/* Dataset pour les systèmes de coordonnées   */
#define ECS_IDEAS_DATASET_SYS_COORD_2420          2420
/* Dataset pour les noeuds */
#define ECS_IDEAS_DATASET_NODES_2411              2411
/* Dataset pour les elements */
#define ECS_IDEAS_DATASET_ELEMENTS_2412           2412
/* Anciens datasets sur les groupes */
#define ECS_IDEAS_DATASET_GROUPS_2430             2430
#define ECS_IDEAS_DATASET_GROUPS_2432             2432
#define ECS_IDEAS_DATASET_GROUPS_2435             2435
#define ECS_IDEAS_DATASET_GROUPS_2452             2452
#define ECS_IDEAS_DATASET_GROUPS_2467             2467
/* Dataset sur les groupes */
#define ECS_IDEAS_DATASET_GROUPS_2477             2477

/* Chaine servant de separateur de dataset */
/*-----------------------------------------*/

#define ECS_IDEAS_SEPARATEUR_DATASET         "    -1\n"

/* Couleurs et groupes */
/*=====================*/

/* Longueur des noms de groupe */
#define ECS_IDEAS_LEN_GROUP_NAME                     40

/* Nombre max d'entites par ligne pour la definition des groupes */
#define ECS_IDEAS_NBR_GRP_ENT_PER_LINE2               2
/* Idem pour les anciens dataset */
#define ECS_IDEAS_NBR_GRP_ENT_PER_LINE4               4

/* Codes des entites appartenant a un groupe */
/*-------------------------------------------*/

#define ECS_IDEAS_TYP_CODE_NODES                      7
#define ECS_IDEAS_TYP_CODE_ELEMENTS                   8

/* Liste des identificateurs de description des elements finis */
/*=============================================================*/

#define ECS_IDEAS_SHELL_LINEAR_TRI                   91  /*  1 */
#define ECS_IDEAS_SHELL_PARABOLIC_TRI                92  /*  2 */
#define ECS_IDEAS_SHELL_CUBIC_TRI                    93  /*  3 */
#define ECS_IDEAS_SHELL_LINEAR_QUAD                  94  /*  4 */
#define ECS_IDEAS_SHELL_PARABOLIC_QUAD               95  /*  5 */
#define ECS_IDEAS_SHELL_CUBIC_QUAD                   96  /*  6 */
#define ECS_IDEAS_SOLID_LINEAR_TETRA                111  /*  7 */
#define ECS_IDEAS_SOLID_PARABOLIC_TETRA             118  /*  8 */
#define ECS_IDEAS_SOLID_LINEAR_WEDGE                112  /*  9 */
#define ECS_IDEAS_SOLID_PARABOLIC_WEDGE             113  /* 10 */
#define ECS_IDEAS_SOLID_CUBIC_WEDGE                 114  /* 11 */
#define ECS_IDEAS_SOLID_LINEAR_BRICK                115  /* 12 */
#define ECS_IDEAS_SOLID_PARABOLIC_BRICK             116  /* 13 */
#define ECS_IDEAS_SOLID_CUBIC_BRICK                 117  /* 14 */

#define ECS_IDEAS_NBR_ELT_TYP                        14

#define ECS_IDEAS_IGNORE                             -1
#define ECS_IDEAS_IGNORE_BEAM                        -2
#define ECS_IDEAS_UNHANDLED                          -3
#define ECS_IDEAS_UNKNOWN                            -4

/* Definition des elements */
/*=========================*/

/* Nombre max de numeros de noeuds par ligne pour la definition des elements */
#define ECS_IDEAS_NBR_NODE_PER_LINE                   8
#define ECS_IDEAS_NBR_ELT_BEAM                        2

/* Tableau donnant la liste des elements `paraboliques' ou `cubiques'
   qui sont transformes en leur equivalent `lineaire'
   ------------------------------------------------------------------*/

#define ECS_IDEAS_ORDER_LINEAR                        1
#define ECS_IDEAS_ORDER_PARABOLIC                     2
#define ECS_IDEAS_ORDER_CUBIC                         3

#define ECS_IDEAS_NBR_MAX_SOM                         8

/*============================================================================
 *  Structures locales
 *============================================================================*/

typedef struct {
  int     label;         /* Identificateur   */
  int     type;          /* 0 = Cartésien    */
  double  transf[3][4];  /* Transformation   */
} _ecs_ideas_sys_coord_t;


typedef struct {
  ecs_int_t     ideas_typ;                      /* Type I-deas de l'element */
  ecs_elt_typ_t ecs_typ;                        /* Type ECS de l'element */
  ecs_int_t     order;                          /* Ordre de l'element */
  ecs_int_t     num_som[ECS_IDEAS_NBR_MAX_SOM]; /* Liste des num de sommet ECS*/
} _ecs_ideas_init_elt_t;

/*============================================================================
 *  Variables globales statiques
 *============================================================================*/

static const int ecs_loc_nbr_max_elt_c = 1000;  /* Nombre initial d'elements  */
static const int ecs_loc_nbr_moy_som_c =    8;  /* Nombre initial de sommets  */
                                                /*  par element en moyenne    */

static const _ecs_ideas_init_elt_t
_ecs_ideas_init_elt_liste_c[ECS_IDEAS_NBR_ELT_TYP] = {

  {                                            /* 1 */
    ECS_IDEAS_SHELL_LINEAR_TRI,
    ECS_ELT_TYP_FAC_TRIA,
    ECS_IDEAS_ORDER_LINEAR,
    { 1, 2, 3 }
  },
  {                                            /* 2 */
    ECS_IDEAS_SHELL_PARABOLIC_TRI,
    ECS_ELT_TYP_FAC_TRIA,
    ECS_IDEAS_ORDER_PARABOLIC,
    { 1, 3, 5 }
  },
  {                                            /* 3 */
    ECS_IDEAS_SHELL_CUBIC_TRI,
    ECS_ELT_TYP_FAC_TRIA,
    ECS_IDEAS_ORDER_CUBIC,
    { 1, 4, 7 }
  },
  {                                            /* 4 */
    ECS_IDEAS_SHELL_LINEAR_QUAD,
    ECS_ELT_TYP_FAC_QUAD,
    ECS_IDEAS_ORDER_LINEAR,
    { 1, 2, 3, 4 }
  },
  {                                            /* 5 */
    ECS_IDEAS_SHELL_PARABOLIC_QUAD,
    ECS_ELT_TYP_FAC_QUAD,
    ECS_IDEAS_ORDER_PARABOLIC,
    { 1, 3, 5, 7 }
  },
  {                                            /* 6 */
    ECS_IDEAS_SHELL_CUBIC_QUAD,
    ECS_ELT_TYP_FAC_QUAD,
    ECS_IDEAS_ORDER_CUBIC,
    { 1, 4, 7, 10 }
  },
  {                                            /* 7 */
    ECS_IDEAS_SOLID_LINEAR_TETRA,
    ECS_ELT_TYP_CEL_TETRA,
    ECS_IDEAS_ORDER_LINEAR,
    { 1, 2, 3, 4 }
  },
  {                                           /*  8 */
    ECS_IDEAS_SOLID_PARABOLIC_TETRA,
    ECS_ELT_TYP_CEL_TETRA,
    ECS_IDEAS_ORDER_PARABOLIC,
    { 1, 3, 5, 10 }
  },
  {                                           /*  9 */
    ECS_IDEAS_SOLID_LINEAR_WEDGE,
    ECS_ELT_TYP_CEL_PRISM,
    ECS_IDEAS_ORDER_LINEAR,
    { 1, 2, 3, 4, 5, 6 }
  },
  {                                           /* 10 */
    ECS_IDEAS_SOLID_PARABOLIC_WEDGE,
    ECS_ELT_TYP_CEL_PRISM,
    ECS_IDEAS_ORDER_PARABOLIC,
    { 1, 3, 5, 10, 12, 14 }
  },
  {                                           /* 11 */
    ECS_IDEAS_SOLID_CUBIC_WEDGE,
    ECS_ELT_TYP_CEL_PRISM,
    ECS_IDEAS_ORDER_CUBIC,
    { 1, 4, 7, 16, 19, 22 }
  },
  {                                           /* 12 */
    ECS_IDEAS_SOLID_LINEAR_BRICK,
    ECS_ELT_TYP_CEL_HEXA,
    ECS_IDEAS_ORDER_LINEAR,
    { 1, 2, 3, 4, 5, 6, 7, 8 },
  },
  {                                           /* 13 */
    ECS_IDEAS_SOLID_PARABOLIC_BRICK,
    ECS_ELT_TYP_CEL_HEXA,
    ECS_IDEAS_ORDER_PARABOLIC,
    { 1, 3, 5, 7, 13, 15, 17, 19 }
  },
  {                                           /* 14 */
    ECS_IDEAS_SOLID_CUBIC_BRICK,
    ECS_ELT_TYP_CEL_HEXA,
    ECS_IDEAS_ORDER_CUBIC,
    { 1, 4, 7, 10, 21, 24, 27, 30 }
  }

};

/*============================================================================
 *  Fonctions privées
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Création d'une table de correspondance des types d'éléments ;
 *  Pour un numéro de type d'élément I-deas donné ityp, la valeur
 *  correspondante à l'indice ityp renvoie l'indice du type correspondant
 *  dans _ecs_ideas_init_elt_liste_c[] si le type est géré, ou une
 *  valeur négative ECS_IDEAS_IGNORE, ECS_IDEAS_UNHANDLED, ou ECS_IDEAS_UNKNOWN
 *  sinon
 *----------------------------------------------------------------------------*/

static ecs_tab_int_t
_ecs_pre_ideas__tab_types(void)
{
  ecs_tab_int_t  tab_types;
  size_t         ind;
  size_t         ityp;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  tab_types.nbr = 232 + 1;
  ECS_MALLOC(tab_types.val, tab_types.nbr, ecs_int_t);

  /* Initialisation */

  for (ityp = 0; ityp < tab_types.nbr; ityp++)
    tab_types.val[ityp] = ECS_IDEAS_UNKNOWN;

  /* Éléments 1 D */

  tab_types.val[ 11] = ECS_IDEAS_IGNORE_BEAM;           /* Rod */
  tab_types.val[ 21] = ECS_IDEAS_IGNORE_BEAM;           /* Linear beam */
  tab_types.val[ 22] = ECS_IDEAS_IGNORE_BEAM;           /* Tapered beam */
  tab_types.val[ 23] = ECS_IDEAS_UNHANDLED;             /* Curved beam */
  tab_types.val[ 24] = ECS_IDEAS_IGNORE_BEAM;           /* Parabolic beam */
  tab_types.val[ 31] = ECS_IDEAS_UNHANDLED;             /* Straight pipe */
  tab_types.val[ 32] = ECS_IDEAS_UNHANDLED;             /* Curved pipe */

  /* Éléments de type "Plane Stress" */

  tab_types.val[ 41] = ECS_IDEAS_SHELL_LINEAR_TRI;
  tab_types.val[ 42] = ECS_IDEAS_SHELL_PARABOLIC_TRI;
  tab_types.val[ 43] = ECS_IDEAS_SHELL_CUBIC_TRI;
  tab_types.val[ 44] = ECS_IDEAS_SHELL_LINEAR_QUAD;
  tab_types.val[ 45] = ECS_IDEAS_SHELL_PARABOLIC_QUAD;
  tab_types.val[ 46] = ECS_IDEAS_SHELL_CUBIC_QUAD;

  /* Éléments de type "Plane Strain" */

  tab_types.val[ 51] = ECS_IDEAS_SHELL_LINEAR_TRI;
  tab_types.val[ 52] = ECS_IDEAS_SHELL_PARABOLIC_TRI;
  tab_types.val[ 53] = ECS_IDEAS_SHELL_CUBIC_TRI;
  tab_types.val[ 54] = ECS_IDEAS_SHELL_LINEAR_QUAD;
  tab_types.val[ 55] = ECS_IDEAS_SHELL_PARABOLIC_QUAD;
  tab_types.val[ 56] = ECS_IDEAS_SHELL_CUBIC_QUAD;

  /* Éléments de type "Plate" */

  tab_types.val[ 61] = ECS_IDEAS_SHELL_LINEAR_TRI;
  tab_types.val[ 62] = ECS_IDEAS_SHELL_PARABOLIC_TRI;
  tab_types.val[ 63] = ECS_IDEAS_SHELL_CUBIC_TRI;
  tab_types.val[ 64] = ECS_IDEAS_SHELL_LINEAR_QUAD;
  tab_types.val[ 65] = ECS_IDEAS_SHELL_PARABOLIC_QUAD;
  tab_types.val[ 66] = ECS_IDEAS_SHELL_CUBIC_QUAD;

  /* Éléments de type "Membrane" */

  tab_types.val[ 71] = ECS_IDEAS_SHELL_LINEAR_TRI;
  tab_types.val[ 72] = ECS_IDEAS_SHELL_PARABOLIC_TRI;
  tab_types.val[ 73] = ECS_IDEAS_SHELL_CUBIC_TRI;
  tab_types.val[ 74] = ECS_IDEAS_SHELL_LINEAR_QUAD;
  tab_types.val[ 75] = ECS_IDEAS_SHELL_PARABOLIC_QUAD;
  tab_types.val[ 76] = ECS_IDEAS_SHELL_CUBIC_QUAD;

  /* Éléments de type "Axisymetric Solid" */

  tab_types.val[ 81] = ECS_IDEAS_SHELL_LINEAR_TRI;
  tab_types.val[ 82] = ECS_IDEAS_SHELL_PARABOLIC_TRI;
  tab_types.val[ 84] = ECS_IDEAS_SHELL_LINEAR_QUAD;
  tab_types.val[ 85] = ECS_IDEAS_SHELL_PARABOLIC_QUAD;

  /* Éléments de type "Thin Shell" */

  tab_types.val[ 91] = ECS_IDEAS_SHELL_LINEAR_TRI;
  tab_types.val[ 92] = ECS_IDEAS_SHELL_PARABOLIC_TRI;
  tab_types.val[ 93] = ECS_IDEAS_SHELL_CUBIC_TRI;
  tab_types.val[ 94] = ECS_IDEAS_SHELL_LINEAR_QUAD;
  tab_types.val[ 95] = ECS_IDEAS_SHELL_PARABOLIC_QUAD;
  tab_types.val[ 96] = ECS_IDEAS_SHELL_CUBIC_QUAD;

  /* Éléments de type "Thick Shell" */

  tab_types.val[101] = ECS_IDEAS_SHELL_LINEAR_TRI;
  tab_types.val[102] = ECS_IDEAS_SHELL_PARABOLIC_TRI;
  tab_types.val[103] = ECS_IDEAS_SHELL_CUBIC_TRI;
  tab_types.val[104] = ECS_IDEAS_SHELL_LINEAR_QUAD;
  tab_types.val[105] = ECS_IDEAS_SHELL_PARABOLIC_QUAD;
  tab_types.val[106] = ECS_IDEAS_SHELL_CUBIC_QUAD;

  /* Éléments de type "Solid" */

  tab_types.val[111] = ECS_IDEAS_SOLID_LINEAR_TETRA;
  tab_types.val[112] = ECS_IDEAS_SOLID_LINEAR_WEDGE;
  tab_types.val[113] = ECS_IDEAS_SOLID_PARABOLIC_WEDGE;
  tab_types.val[114] = ECS_IDEAS_SOLID_CUBIC_WEDGE;
  tab_types.val[115] = ECS_IDEAS_SOLID_LINEAR_BRICK;
  tab_types.val[116] = ECS_IDEAS_SOLID_PARABOLIC_BRICK;
  tab_types.val[117] = ECS_IDEAS_SOLID_CUBIC_BRICK;
  tab_types.val[118] = ECS_IDEAS_SOLID_PARABOLIC_TETRA;

  /* Autres éléments, ignorés */

  tab_types.val[121] = ECS_IDEAS_IGNORE;
  tab_types.val[122] = ECS_IDEAS_IGNORE;
  tab_types.val[136] = ECS_IDEAS_IGNORE;
  tab_types.val[137] = ECS_IDEAS_IGNORE;
  tab_types.val[138] = ECS_IDEAS_IGNORE;
  tab_types.val[139] = ECS_IDEAS_IGNORE;
  tab_types.val[141] = ECS_IDEAS_IGNORE;
  tab_types.val[151] = ECS_IDEAS_IGNORE;
  tab_types.val[152] = ECS_IDEAS_IGNORE;
  tab_types.val[161] = ECS_IDEAS_IGNORE;  /* Lumped Mass */
  tab_types.val[171] = ECS_IDEAS_IGNORE;  /* Axisymetric Linear Shell */
  tab_types.val[172] = ECS_IDEAS_IGNORE;  /* Axisymetric Parabolic Shell */
  tab_types.val[181] = ECS_IDEAS_IGNORE;  /* Constraint */
  tab_types.val[191] = ECS_IDEAS_IGNORE;  /* Plastic Cold Runner */
  tab_types.val[192] = ECS_IDEAS_IGNORE;  /* Plastic Hot Runner */
  tab_types.val[193] = ECS_IDEAS_IGNORE;  /* Plastic Water Line */
  tab_types.val[194] = ECS_IDEAS_IGNORE;  /* Plastic Fountain */
  tab_types.val[195] = ECS_IDEAS_IGNORE;  /* Plastic Baffle */
  tab_types.val[196] = ECS_IDEAS_IGNORE;  /* Plastic Rod Heater */
  tab_types.val[201] = ECS_IDEAS_IGNORE;  /* Linear node-to-node interface */
  tab_types.val[202] = ECS_IDEAS_IGNORE;  /* Linear edge-to-edge interface */
  tab_types.val[203] = ECS_IDEAS_IGNORE;  /* Parab. edge-to-edge interface */
  tab_types.val[204] = ECS_IDEAS_IGNORE;  /* Linear face-to-face interface */
  tab_types.val[208] = ECS_IDEAS_IGNORE;  /* Parab. face-to-face interface */
  tab_types.val[212] = ECS_IDEAS_IGNORE;  /* Linear axisymetric interface */
  tab_types.val[213] = ECS_IDEAS_IGNORE;  /* Parab. axisymetric interface */
  tab_types.val[221] = ECS_IDEAS_IGNORE;  /* Linear rigid surface */
  tab_types.val[222] = ECS_IDEAS_IGNORE;  /* Parabolic rigid surface */
  tab_types.val[231] = ECS_IDEAS_IGNORE;  /* Axisym. linear rigid surface */
  tab_types.val[232] = ECS_IDEAS_IGNORE;  /* Axisym. parab. rigid surface */

  /* Renumérotation des types gérés */

  for (ityp = 0; ityp < tab_types.nbr; ityp++) {

    if (tab_types.val[ityp] > 0) {

      for (ind = 0;
           (   ind < ECS_IDEAS_NBR_ELT_TYP
            && (  _ecs_ideas_init_elt_liste_c[ind].ideas_typ
                != tab_types.val[ityp]));
           ind++);

      assert(ind < ECS_IDEAS_NBR_ELT_TYP);
      tab_types.val[ityp] = ind;

    }

  }

  /* Renvoi du tableau crée et initialisé */

  return tab_types;
}

/*----------------------------------------------------------------------------
 *  Fonction utilitaire pour la lecture des réels d'un fichier I-DEAS
 *
 *  La chaine contenant un réel avec (ou non) un exposant `d' ou `D'
 *  est convertie en réel avec un exposant `e'
 *----------------------------------------------------------------------------*/

static ecs_coord_t
_ecs_pre_ideas__transf_expo(char *chaine)
{
  char        *s;
  char         s_atof[ECS_LOC_LNG_MAX_CHAINE_IDEAS];
  ecs_coord_t  f_atof;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  s = s_atof;

  /* Copie de la partie de la chaine en entree `chaine'              */
  /*  qui precede le caractere `D' ou `d'                            */
  /*  ou copie de toute la chaine si `D' ou `d' ne sont pas presents */

  while(*chaine != 'D' && *chaine != 'd' && *chaine != '\0') {
    *s = *chaine;
    s++;
    chaine++;
  }

  /* Toute la chaine en entree a-t-elle ete parcourue ? */

  if (*chaine != '\0') {

    /* La chaine en entree n'a pas ete entierement parcourue :     */
    /*  le caractere `D' ou `d' est remplace par `e' dans sa copie */

    *s = 'e';

    /* Le reste de la chaine est copie tel quel */

    while ((*++s = *++chaine));

  }
  else {

    /* Toute la chaine a ete parcourue.      */
    /* On ajoute le `\0' final dans la copie */

    *s = *chaine;

  }

  /* L'exposant `D' ou `d' a ete remplace par `e' :                  */
  /*  on peut maintenant extraire le reel de la chaine de caracteres */

  sscanf(s_atof, " %lf", &f_atof);

  return f_atof;
}


/*----------------------------------------------------------------------------
 *  Lecture des systèmes de coordonnées
 *----------------------------------------------------------------------------*/

static void
_ecs_pre_ideas__lit_sys_coord(ecs_file_t               *fic_maillage,
                              int                      *num_ligne,
                              ecs_tab_int_t            *ind_coo_sys,
                              _ecs_ideas_sys_coord_t  **coo_sys_ref)
{
  ecs_int_t    retour;                /* Retour fonctions `scan' pour test    */
  ecs_int_t    num_rec;               /* Numéro local de l'enregistrement     */
  ecs_int_t    iloc;
  ecs_int_t    jloc;

  bool         debut_part;
  bool         erreur_decod;
  bool         sys_trivial;

  _ecs_ideas_sys_coord_t      coo_sys_loc;
  ecs_int_t                   nbr_part;
  ecs_int_t                   nbr_sys_triv;
  ecs_int_t                   nbr_sys_cyl;
  ecs_int_t                   nbr_sys_sphr;
  size_t                      num_label_max;

  char         chaine[ECS_LOC_LNG_MAX_CHAINE_IDEAS];          /* Ligne lue   */

  /* Variables I-DEAS lues */
  int          ideas_part_uid;           /* Numéro du "part" I-DEAS           */
  int          ideas_sys_coo_num;        /* Numéro du système de coordonnées  */
  int          ideas_sys_coo_type;       /* Type du système de coordonnées    */
  int          ideas_sys_coo_color;      /* Couleur du système                */

  char         ideas_ch_val[3][ECS_LOC_LNG_MAX_CHAINE_IDEAS];
                                         /* Chaines receptrices des valeurs   */
                                         /* avant transformation en réel      */

  ecs_coord_t  transf_triv[3][4] = {{1.0, 0.0, 0.0, 0.0},
                                    {0.0, 1.0, 0.0, 0.0},
                                    {0.0, 0.0, 1.0, 0.0}};

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


  /* Initialisations */
  /*=================*/

  nbr_part = 0;
  nbr_sys_triv = 0;
  nbr_sys_cyl  = 0;
  nbr_sys_sphr = 0;
  num_label_max = 0;

  ind_coo_sys->nbr = 5;
  ECS_MALLOC(ind_coo_sys->val, ind_coo_sys->nbr, ecs_int_t);

  *coo_sys_ref = NULL;

  num_rec = 1;

  debut_part = false;
  erreur_decod = false;


  /*========================================================*/
  /* Tant qu'on n'a pas atteint la fin de fichier           */
  /* et que la ligne lue n'est pas le separateur de dataset */
  /*========================================================*/

  while (ecs_file_gets_try(chaine, ECS_LOC_LNG_MAX_CHAINE_IDEAS,
                           fic_maillage, num_ligne) != NULL
         && strcmp(chaine, ECS_IDEAS_SEPARATEUR_DATASET) != 0) {


    /* Lecture d'un système de coordonnées */
    /*=====================================*/

    /*
     * On ne sait pas à priori si on change de "part" (enregistrements
     * 1 et 2 non répétés pour systèmes d'un même "part") ; dans un cas
     * on lit les enregistrements 1 à 8, dans l'autre on ne lit que les
     * enregistrements 3 à 8.
     */

    switch (num_rec) {

    case 1: /* UID de part (format 1I10) */
    case 3: /* label, type, et couleur du système (format 3I10) */

      /* On vérifie si on n'est pas déjà à l'enregistrement 3 :          */

      retour = sscanf(chaine,
                      " %d %d %d",
                      &ideas_sys_coo_num,
                      &ideas_sys_coo_type,
                      &ideas_sys_coo_color);

      if (retour == 3) {

        coo_sys_loc.label = ideas_sys_coo_num;
        coo_sys_loc.type  = ideas_sys_coo_type;

        if (num_rec == 1) /* On a bien lu l'enregistrement 3, et non 1 */
          num_rec = 3;

      }
      else if (num_rec == 1) {

        /* UID de part (format 3I10) */
        retour = sscanf(chaine, " %d", &ideas_part_uid);

        if (retour != 1)
          erreur_decod = true;
        else
          debut_part = true;

      }

      else
        erreur_decod = true;

      break;

    case 2: /* Nom de part associé           (format 40A2) inutile ici */
    case 4: /* Nom de système de coordonnées (format 40A2) inutile ici */

      break;

    case 5: /* 1ère rangée de la matrice de transformation (format 1P3D25.16) */
    case 6: /* 2ème rangée de la matrice de transformation (format 1P3D25.16) */
    case 7: /* 3ème rangée de la matrice de transformation (format 1P3D25.16) */
    case 8: /* 4ème rangée de la matrice de transformation (format 1P3D25.16) */

      /* I-DEAS peut sortir des coordonnées avec un exposant D        */
      /* qui ne peuvent pas etre relues de maniere portable par le C. */
      /* On lit donc les coordonnées en les placant dans des chaînes  */

    /* Format Ideas : 1P3D25.16 */

    retour = sscanf(chaine, " %s %s %s",
                    ideas_ch_val[0], ideas_ch_val[1], ideas_ch_val[2]);

    if (retour != 3)
      erreur_decod = true;

    else {

    /* On transforme les coordonnées sous forme de chaîne          */
    /* en substituant l'exposant D par un exposant E et on renvoie */
    /* le réel correspondant                                       */

    for (jloc = 0; jloc < 3; jloc++)
      coo_sys_loc.transf[jloc][num_rec - 5]
        = _ecs_pre_ideas__transf_expo(ideas_ch_val[jloc]);

    }

    break;

    default:
      assert(num_rec > 0 && num_rec < 9);

    }

    if (erreur_decod == true)
      ecs_error(__FILE__, __LINE__, errno,
                _("Error reading line %d of file \"%s\"."),
                *num_ligne, ecs_file_get_name(fic_maillage));

    num_rec += 1;


    /* On a terminé pour le système de coordonnées en cours */
    /*------------------------------------------------------*/

    if (num_rec > 8) {

      num_rec = 1;

      /* En fait, les coordonnées des points sont décrites en fonction
         du système de coordonnées principal attaché à la part,
         et non du système de coordonnées local */

      if (debut_part == true) {

        nbr_part += 1;
        ECS_REALLOC(*coo_sys_ref, nbr_part, _ecs_ideas_sys_coord_t);

        (*coo_sys_ref)[nbr_part - 1] = coo_sys_loc;

        debut_part = false;

        /* Repérage des systèmes de coordonnées "triviaux"
           ou non cartésiens */

        sys_trivial = true;

        if (coo_sys_loc.type == 0) {
          for (jloc = 0; jloc < 3; jloc++) {
            for (iloc = 0; iloc < 4; iloc++) {
              if (ECS_ABS(  coo_sys_loc.transf[jloc][iloc]
                          - transf_triv[jloc][iloc])
                  > ECS_REAL_PRECISION)
                sys_trivial = false;
            }
          }
        }
        else {
          sys_trivial = false;
          if (coo_sys_loc.type == 1)
            nbr_sys_cyl += 1;
          else if (coo_sys_loc.type == 2)
            nbr_sys_sphr += 1;
        }

        if (sys_trivial == true)
          nbr_sys_triv += 1;

      }

      num_label_max = ECS_MAX(coo_sys_loc.label, (ecs_int_t)num_label_max);

      if (num_label_max > ind_coo_sys->nbr) {
        ind_coo_sys->nbr = num_label_max*2;
        ECS_REALLOC(ind_coo_sys->val, ind_coo_sys->nbr, ecs_int_t);
      }

      ind_coo_sys->val[coo_sys_loc.label - 1] = nbr_part - 1;

    }

  }

  ind_coo_sys->nbr = num_label_max;
  ECS_REALLOC(ind_coo_sys->val, ind_coo_sys->nbr, ecs_int_t);

  /* Les systèmes de coordonnées ne pris en compte que si nécessaire */

  if (nbr_sys_cyl > 0 || nbr_sys_sphr >0) {
    printf(_("    %d cylindrical system(s) and "
             "%d spherical system(s) not converted;\n"),
           (int)nbr_sys_cyl, (int)nbr_sys_sphr);

    printf(_("    probable problems depending on I-deas "
             "construction and export\n\n"));
  }

  if (getenv("CS_PREPROCESS_IGNORE_IDEAS_COO_SYS") != NULL) {
    if (atoi(getenv("CS_PREPROCESS_IGNORE_IDEAS_COO_SYS")) > 0) {

      if (nbr_part - nbr_sys_triv > 0) {

        printf(_("  Ignored coordinate systems:\n"));
        printf(_("    %d non-identity transformation(s) ignored:\n"),
               (int)(nbr_part - nbr_sys_triv));

        nbr_sys_triv = nbr_part;
      }
    }
  }

  if (nbr_sys_triv == nbr_part) {
    ind_coo_sys->nbr = 0;
    ECS_FREE(ind_coo_sys->val);
    ECS_FREE(*coo_sys_ref);
  }

}


/*----------------------------------------------------------------------------
 *  Lecture des coordonnées des noeuds
 *----------------------------------------------------------------------------*/

static void
_ecs_pre_ideas__lit_nodes(ecs_maillage_t            *maillage,
                          ecs_file_t                *fic_maillage,
                          int                       *num_ligne,
                          ecs_int_t                **som_val_label,
                          bool                      *bool_label_som_a_trier,
                          ecs_tab_int_t              ind_coo_sys,
                          _ecs_ideas_sys_coord_t     coo_sys_ref[])
{
  ecs_int_t    retour;                /* Retour des fnctions `scan' pour test */
  char       chaine[ECS_LOC_LNG_MAX_CHAINE_IDEAS];             /* Ligne lue   */

  /* Variables Ideas lues */
  int          ideas_num_nod;            /* Numero  Ideas du noeud            */
  int          ideas_sys_coo;            /* Sys. coo. de référence du noeud   */

  char       ideas_ch_coord[3][ECS_LOC_LNG_MAX_CHAINE_IDEAS];
                                         /* Chaines receptrices des coords    */
                                         /* avant transformation en reel      */

  ecs_coord_t  coord[3];                 /* Coordonnees reelles du sommet     */
  ecs_int_t    max_som;                  /* Dim d'allocation de tab. locaux   */
  ecs_int_t    nbr_som_max;              /* Nombre initial de sommets         */

  ecs_int_t    icoo;                     /* Indice de boucle sur les comp.    */
                                         /* des coordonnees (X, Y et Z)       */

  /* Stockage avant transfert */
  /*--------------------------*/

  ecs_int_t    cpt_som = 0;              /* Compteur des sommets lus          */

  ecs_coord_t * som_val_coord ;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


  /* Initialisations */
  /*=================*/

  nbr_som_max = ecs_loc_nbr_max_elt_c * ecs_loc_nbr_moy_som_c;

  max_som   = nbr_som_max;

  ECS_MALLOC(som_val_coord,   max_som * 3, ecs_coord_t);
  ECS_MALLOC(*som_val_label,  max_som,   ecs_int_t);


  /*========================================================*/
  /* Tant qu'on n'a pas atteint la fin de fichier           */
  /* et que la ligne lue n'est pas le separateur de dataset */
  /*========================================================*/

  while (ecs_file_gets_try(chaine, ECS_LOC_LNG_MAX_CHAINE_IDEAS,
                           fic_maillage, num_ligne) != NULL
         && strcmp(chaine, ECS_IDEAS_SEPARATEUR_DATASET) != 0) {


    /* Lecture des caracteristiques du noeud */
    /*=======================================*/

    /* Format Ideas : 4I10 */

    retour = sscanf(chaine,
                    " %d %d",
                    &ideas_num_nod, &ideas_sys_coo);

    if (retour != 2)
        ecs_error(__FILE__, __LINE__, errno,
                  _("Error reading line %d of file \"%s\"."),
                  *num_ligne, ecs_file_get_name(fic_maillage));

    /* Lecture des coordonnees du sommet */
    /*===================================*/

    ecs_file_gets(chaine, ECS_LOC_LNG_MAX_CHAINE_IDEAS,
                  fic_maillage, num_ligne);

    /* Ideas peut sortir des coordonnees avec un exposant D         */
    /* qui ne peuvent pas etre relues de maniere portable par le C. */
    /* On lit donc les coordonnees en les placant dans des chaines  */

    /* Format Ideas : 1P3D25.16 */

    retour = sscanf(chaine, " %s %s %s",
                    ideas_ch_coord[0], ideas_ch_coord[1], ideas_ch_coord[2]);

    if (retour != 3)
      ecs_error(__FILE__, __LINE__, errno,
                _("Error reading line %d of file \"%s\"."),
                *num_ligne, ecs_file_get_name(fic_maillage));

    /* On transforme les coordonnees sous forme de chaine          */
    /* en substituant l'exposant D par un exposant E et on renvoit */
    /* le reel correspondant                                       */

    for (icoo = 0; icoo < 3; icoo++) {
      coord[icoo] = _ecs_pre_ideas__transf_expo(ideas_ch_coord[icoo]);
    }


    /* Stockage des valeurs lues avant transfert dans la structure `maillage' */
    /*========================================================================*/

    /* Reallocation eventuelles des tableaux locaux */
    /*----------------------------------------------*/

    if (max_som   <= cpt_som) {

      max_som = ECS_MAX(max_som * 2, cpt_som + 1);

      ECS_REALLOC(som_val_coord,   max_som * 3, ecs_coord_t);
      ECS_REALLOC(*som_val_label,  max_som, ecs_int_t);

    }


    /* Coordonnees du sommet lu */
    /*--------------------------*/

    if (ind_coo_sys.nbr > 0) {

      ideas_sys_coo = ind_coo_sys.val[ideas_sys_coo - 1];

      for (icoo = 0; icoo < 3; icoo++) {
        som_val_coord[cpt_som * 3 + icoo]
          = (  coo_sys_ref[ideas_sys_coo].transf[icoo][0] * coord[0]
             + coo_sys_ref[ideas_sys_coo].transf[icoo][1] * coord[1]
             + coo_sys_ref[ideas_sys_coo].transf[icoo][2] * coord[2]
             + coo_sys_ref[ideas_sys_coo].transf[icoo][3]);
      }

    }
    else {

      for (icoo = 0; icoo < 3; icoo++)
        som_val_coord[cpt_som * 3 + icoo] = coord[icoo];

    }

    /* Etiquette   du sommet lu */
    /*--------------------------*/

    if ((*bool_label_som_a_trier == false) &&
        (cpt_som > 0) && ((*som_val_label)[cpt_som - 1] > ideas_num_nod))
      *bool_label_som_a_trier = true;

    (*som_val_label)[cpt_som] = ideas_num_nod;

    /* Incrementation du compteur de sommets lus */
    /*===========================================*/

    cpt_som++;

  } /* Fin : tant qu'il y a une ligne a lire et                       */
    /*       tant que la ligne lue n'est pas le separateur de dataset */

  if (strcmp(chaine, ECS_IDEAS_SEPARATEUR_DATASET) != 0)
    ecs_file_read_check_error(fic_maillage, *num_ligne);

  /* else : on a lu jusqu'au separateur de dataset */

  /* Reallocations des tableaux locaux */
  /*===================================*/

  ECS_REALLOC(som_val_coord,   cpt_som * 3, ecs_coord_t);
  ECS_REALLOC(*som_val_label,  cpt_som,     ecs_int_t);

  /* Transfert des valeurs lues dans la structure d'entite de maillage */
  /*===================================================================*/

  ecs_maillage_pre__cree_som(maillage, cpt_som, som_val_coord);
}


/*----------------------------------------------------------------------------
 *  Lecture de la table de connectivité
 *----------------------------------------------------------------------------*/

static void
_ecs_pre_ideas__lit_elements(ecs_maillage_t   *maillage,
                             ecs_file_t       *fic_maillage,
                             int              *num_ligne,
                             int               nbr_som,
                             ecs_int_t       **som_val_label,
                             ecs_tab_int_t     tab_label_ent[],
                             bool              bool_label_ent_a_trier[])
{
  int        retour;                  /* Retour fonctions `scan' pour test    */
  char       chaine[ECS_LOC_LNG_MAX_CHAINE_IDEAS];             /* Ligne lue   */
  char      *ssch;                    /* Sous-chaine de lecture               */

  /* Variables Ideas lues */
  int          ideas_num_elt;         /* Numero         Ideas de l'element    */
  int          ideas_fe_des_id;       /* Identificateur Ideas de l'element    */
  int          ideas_color;           /* Couleur        Ideas de l'element    */
  int          ideas_nbr_nod_elt;     /* Nombre de noeuds     de l'element    */
  int          ideas_nod_elt[ECS_IDEAS_NBR_MAX_SOM];
                                      /* Numeros des noeuds de l'element      */

  int          nbr_som_elt;           /* Nb de noeuds a lire = nbr de sommets */
  int          nbr_nod_elt_reste;     /* Nb de noeuds de l'elt restant a lire */
  int          nbr_nod_ligne;         /* Nb de noeuds par ligne lue           */
  ecs_int_t    ent_num;               /* Numero de l'entite concernee         */
  size_t       max_som_ent[ECS_N_ENTMAIL];   /* Longueur alloc tableaux loc.  */
  size_t       max_elt_ent[ECS_N_ENTMAIL];   /* Longueur alloc tableaux loc.  */
  size_t       num_pos;               /* Nouvel indice de stockage            */

  ecs_int_t    icoul;                 /* Indice de boucle sur les couleurs    */
  ecs_int_t    ient;                  /* Indice de boucle sur les entites     */
  ecs_int_t    inod;                  /* Indice de boucle sur les noeuds      */
  ecs_int_t    inod_lig;              /* Indice de boucle des noeuds sur      */
                                      /*  1 ligne de lecture                  */
  ecs_int_t    inum_som;              /* Indice Ideas des sommets             */
  ecs_int_t    isom;                  /* Indice de boucle sur les sommets     */
  ecs_int_t    ityp;                  /* Boucle sur les types d'elements      */
  ecs_int_t    cpt_coul_ent[ECS_N_ENTMAIL];   /* Compteur de couleurs         */
  ecs_int_t   *val_coul_ent[ECS_N_ENTMAIL];   /* Tableau valeurs des couleurs */
  ecs_size_t  *cpt_elt_coul_ent[ECS_N_ENTMAIL];

  ecs_tab_int_t  corresp_typ;         /* Tableau de correspondance des types  */

  /* Stockage avant transfert */
  /*--------------------------*/

  size_t       cpt_elt_ent        [ECS_N_ENTMAIL]; /* Nombre d'elems/entite */

  ecs_size_t  *elt_pos_som_ent    [ECS_N_ENTMAIL]; /* Positions numeros som */
  ecs_int_t   *elt_val_som_ent    [ECS_N_ENTMAIL]; /* Numeros des sommets   */
  ecs_int_t   *elt_val_label_ent  [ECS_N_ENTMAIL]; /* Etiquettes            */
  ecs_int_t   *elt_val_color_ent  [ECS_N_ENTMAIL]; /* Couleurs              */

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /*====================================================*/
  /* Initialisations et allocations des tableaux locaux */
  /*====================================================*/

  /* Attention au decalage de `1' !!!         */
  /* On n'alloue pas les tableaux locaux pour */
  /* `ECS_ENTMAIL_DEB = ECS_ENTMAIL_SOM'      */

  for (ient = ECS_ENTMAIL_FAC; ient < ECS_N_ENTMAIL; ient++) {

    cpt_elt_ent        [ient] = 0;

    cpt_coul_ent       [ient] = 0;
    cpt_elt_coul_ent   [ient] = NULL;
    val_coul_ent       [ient] = NULL;

    max_som_ent        [ient] = ecs_loc_nbr_max_elt_c * ecs_loc_nbr_moy_som_c;
    max_elt_ent        [ient] = ecs_loc_nbr_max_elt_c;

  }

  corresp_typ = _ecs_pre_ideas__tab_types();

  /*========================================================*/
  /* Tant qu'on n'a pas atteint la fin de fichier           */
  /* et que la ligne lue n'est pas le separateur de dataset */
  /*========================================================*/

#define ECS_FCT_TYP(ityp) _ecs_ideas_init_elt_liste_c[ityp].ecs_typ

  while (ecs_file_gets_try(chaine, ECS_LOC_LNG_MAX_CHAINE_IDEAS,
                           fic_maillage, num_ligne) != NULL
         && strcmp(chaine, ECS_IDEAS_SEPARATEUR_DATASET) != 0){

    /* Lecture des caracteristiques de l'element */
    /*===========================================*/

    /* Format Ideas : 6I10 */

    retour = sscanf(chaine,
                    " %d %d %*d %*d %d %d",
                    &ideas_num_elt,
                    &ideas_fe_des_id,
                    &ideas_color,
                    &ideas_nbr_nod_elt);

    if (retour != 4)
      ecs_error(__FILE__, __LINE__, errno,
                _("Error reading line %d of file \"%s\"."),
                *num_ligne, ecs_file_get_name(fic_maillage));

    /* Est-ce un element dont le type Ideas est reconnu ?              */
    /* Est-ce un element `non lineaire' (`cubique' ou `parabolique') ? */
    /*-----------------------------------------------------------------*/

    if (ideas_fe_des_id < (int)(corresp_typ.nbr))
      ityp = corresp_typ.val[ideas_fe_des_id];
    else
      ityp = ECS_IDEAS_UNKNOWN;

    if (ityp == ECS_IDEAS_UNKNOWN || ityp == ECS_IDEAS_UNHANDLED) {

      /* Le type Ideas de l'element n'est pas reconnu */

      if (ityp == ECS_IDEAS_UNKNOWN)
        ecs_error(__FILE__, __LINE__, 0,
                  _("Error reading an I-deas universal file:\n"
                    "at line %d of file \"%s\".\n"
                    "Type identifier <%d> of element <%d> not recognized."),
                  *num_ligne, ecs_file_get_name(fic_maillage),
                  (int)ideas_fe_des_id, (int)ideas_num_elt);

      else /* if (ityp == ECS_IDEAS_UNHANDLED) */
        ecs_error(__FILE__, __LINE__, 0,
                  _("Error reading an I-deas universal file:\n"
                    "at line %d of file \"%s\".\n"
                    "Type identifier <%d> of element <%d> recognized "
                    "but not handled."),
                  *num_ligne, ecs_file_get_name(fic_maillage),
                  (int)ideas_fe_des_id, (int)ideas_num_elt);

    }
    /* else :  le type Ideas de l'element est reconnu */

    /* Est-ce un element `poutre' (`beam') ? */
    /*---------------------------------------*/

    if (ityp == ECS_IDEAS_IGNORE_BEAM) {

      /* On saute la ligne qui contient les tables propres aux `beam' */

      /* Format Ideas : 3I10 */

      ecs_file_gets(chaine, ECS_LOC_LNG_MAX_CHAINE_IDEAS,
                    fic_maillage, num_ligne);
    }

    /* Si ce type d'element est à ignorer, on passe au suivant */

    if (ityp < 0) {
      ecs_file_gets(chaine, ECS_LOC_LNG_MAX_CHAINE_IDEAS,
                    fic_maillage, num_ligne);
      continue;
    }

    /* Lecture des numeros de noeuds constituant l'element */
    /*=====================================================*/

    /* Format Ideas : 8I10 (pour chaque ligne) */

    /* Nombre de noeuds a lire =              */
    /*  nombre de noeuds qui sont des sommets */
    nbr_som_elt = ecs_fic_elt_typ_liste_c[ECS_FCT_TYP(ityp)].nbr_som;

    nbr_nod_elt_reste = ideas_nbr_nod_elt;

    isom     = 0;
    inod     = 0;
    inum_som = 0;

    do {

      ecs_file_gets(chaine, ECS_LOC_LNG_MAX_CHAINE_IDEAS,
                    fic_maillage, num_ligne);

      /* Lecture des numeros de noeuds pour la ligne courante */

      nbr_nod_ligne
        = ECS_MIN(nbr_nod_elt_reste, ECS_IDEAS_NBR_NODE_PER_LINE);

      ssch = strtok(chaine, " ");

      for (inod_lig = 0; inod_lig <  nbr_nod_ligne; inod_lig++) {

        if (inod == (_ecs_ideas_init_elt_liste_c[ityp].num_som[inum_som]-1)){
          retour = sscanf(ssch,"%d", &ideas_nod_elt[isom++]);
          ssch   = strtok(NULL, " ");
          inum_som++;
        }
        else {
          retour = sscanf(ssch, "%*d");
          ssch   = strtok(NULL, " ");
        }

        inod++;

      }
      if (ssch != NULL)
        ecs_error(__FILE__, __LINE__, errno,
                  _("Error reading line %d of file \"%s\"."),
                  *num_ligne, ecs_file_get_name(fic_maillage));

      nbr_nod_elt_reste = nbr_nod_elt_reste - nbr_nod_ligne;

    } while (nbr_nod_elt_reste != 0);

    if (isom != nbr_som_elt)
      ecs_error(__FILE__, __LINE__, errno,
                _("Error reading line %d of file \"%s\"."),
                *num_ligne, ecs_file_get_name(fic_maillage));

    /* Stockage des valeurs lues avant transfert dans la structure `maillage' */
    /*========================================================================*/

    /* Identification de l'entite concernee */
    /*--------------------------------------*/

    ent_num =
      ecs_maillage_pre__ret_typ_geo(_ecs_ideas_init_elt_liste_c[ityp].ecs_typ);

    if (cpt_elt_ent[ent_num] != 0) {

      /* Reallocations eventuelles des tableaux locaux */
      /*-----------------------------------------------*/

      if (max_elt_ent[ent_num] <= cpt_elt_ent[ent_num]) {

        max_elt_ent[ent_num] = ECS_MAX(max_elt_ent[ent_num] * 2,
                                       cpt_elt_ent[ent_num] + 1);

        ECS_REALLOC(elt_val_label_ent[ent_num],   max_elt_ent[ent_num],
                    ecs_int_t);

        ECS_REALLOC(elt_val_color_ent[ent_num],   max_elt_ent[ent_num],
                    ecs_int_t);

        ECS_REALLOC(elt_pos_som_ent[ent_num],     max_elt_ent[ent_num] + 1,
                    ecs_size_t);

      }

    }
    else { /* cpt_elt_ent[ent_num] == 0 */

      /* C'est la premiere fois que l'entite est concernee */
      /* On initialise les allocations des tableaux        */
      /*---------------------------------------------------*/

      ECS_MALLOC(elt_pos_som_ent[ent_num],
                 max_elt_ent[ent_num] + 1,     ecs_size_t);
      ECS_MALLOC(elt_val_som_ent[ent_num],
                 max_som_ent[ent_num],         ecs_int_t);
      ECS_MALLOC(elt_val_label_ent[ent_num],
                 max_elt_ent[ent_num],         ecs_int_t);
      ECS_MALLOC(elt_val_color_ent[ent_num],
                 max_elt_ent[ent_num],         ecs_int_t);

      elt_pos_som_ent[ent_num][0] = 1;

    }

    num_pos = elt_pos_som_ent[ent_num][cpt_elt_ent[ent_num]] - 1 + nbr_som_elt;

    if (max_som_ent[ent_num] <= num_pos) {
      max_som_ent[ent_num] = ECS_MAX(max_som_ent[ent_num] * 2, num_pos + 1);
      ECS_REALLOC(elt_val_som_ent[ent_num], max_som_ent[ent_num], ecs_int_t);
    }

    /* Determination de la position des numeros de sommets */
    /* du prochain element                                 */
    /*-----------------------------------------------------*/

    elt_pos_som_ent[ent_num][cpt_elt_ent[ent_num] + 1] =
      elt_pos_som_ent[ent_num][cpt_elt_ent[ent_num]] + nbr_som_elt;

    /* Connectivite de l'element par ses numeros de sommets */
    /*------------------------------------------------------*/

    for (inod = 0; inod <  nbr_som_elt; inod++) {

      elt_val_som_ent
        [ent_num]
        [elt_pos_som_ent[ent_num][cpt_elt_ent[ent_num]] - 1 + inod]
        = ideas_nod_elt[inod];

    }

    /* Etiquette de l'element lu */
    /*---------------------------*/

    if (((bool_label_ent_a_trier[ent_num]) == false) &&
        (cpt_elt_ent[ent_num] > 0) &&
        (elt_val_label_ent[ent_num][cpt_elt_ent[ent_num] - 1] > ideas_num_elt))
        bool_label_ent_a_trier[ent_num] = true;

    elt_val_label_ent[ent_num][cpt_elt_ent[ent_num]] = ideas_num_elt;

    /* Couleur   de l'element lu */
    /*---------------------------*/

    icoul = 0;
    while (icoul < cpt_coul_ent[ent_num]               &&
           val_coul_ent[ent_num][icoul] != ideas_color   )
      icoul++;

    if (icoul == cpt_coul_ent[ent_num]) {

      /* La valeur de la couleur n'a pas encore ete stockee */

      ECS_REALLOC(val_coul_ent[ent_num],     cpt_coul_ent[ent_num] + 1,
                  ecs_int_t);
      ECS_REALLOC(cpt_elt_coul_ent[ent_num], cpt_coul_ent[ent_num] + 1,
                  ecs_size_t);
      cpt_elt_coul_ent[ent_num][icoul] = 0;
      val_coul_ent[ent_num][icoul] = ideas_color;
      cpt_coul_ent[ent_num]++;

    }

    cpt_elt_coul_ent[ent_num][icoul]++;
    elt_val_color_ent[ent_num][cpt_elt_ent[ent_num]] = icoul + 1;

    /* Incrementation du nombre d'elements lus */
    /*=========================================*/

    cpt_elt_ent[ent_num]++;


  } /* Fin : tant qu'il y a une ligne a lire et
     *       tant que la ligne lue n'est pas le separateur de dataset */

#undef ECS_FCT_TYP


  if (strcmp(chaine, ECS_IDEAS_SEPARATEUR_DATASET) != 0)

    ecs_file_read_check_error(fic_maillage,
                              *num_ligne);

  /* else : on a lu jusqu'au separateur de dataset */

  corresp_typ.nbr = 0;
  ECS_FREE(corresp_typ.val);

  /* Reallocations des tableaux locaux et mise à jour des références */
  /*=================================================================*/

  for (ient = ECS_ENTMAIL_FAC; ient < ECS_N_ENTMAIL; ient++) {

    if (cpt_elt_ent[ient] != 0) {
      ECS_REALLOC(elt_pos_som_ent[ient],
                  cpt_elt_ent[ient] + 1,     ecs_size_t);
      ECS_REALLOC(elt_val_som_ent[ient],
                  elt_pos_som_ent[ient][cpt_elt_ent[ient]] - 1, ecs_int_t);
      ECS_REALLOC(elt_val_label_ent[ient],
                  cpt_elt_ent[ient],         ecs_int_t);
      ECS_REALLOC(elt_val_color_ent[ient],
                  cpt_elt_ent[ient],         ecs_int_t);

      ecs_maillage_pre__label_en_indice
        (nbr_som,
         elt_pos_som_ent[ient][cpt_elt_ent[ient]] - 1,
         *som_val_label,
         elt_val_som_ent[ient]);
    }
  }

  ECS_FREE(*som_val_label);

  /* On conserve les listes des labels pour la lecture des groupes */
  /*===============================================================*/

  for (ient = ECS_ENTMAIL_FAC; ient < ECS_N_ENTMAIL; ient++) {

    /* S'il y a au moins un element de cette entite de maillage */

    if (cpt_elt_ent[ient] != 0) {

      tab_label_ent[ient].nbr = cpt_elt_ent[ient];

      tab_label_ent[ient].val = elt_val_label_ent[ient];
      elt_val_label_ent[ient] = NULL;

    }
  }

  /* Transfert des valeurs lues dans les structures d'entite de maillage */
  /*=====================================================================*/

  ecs_maillage_pre__cree_elt(maillage,
                             cpt_elt_ent,
                             elt_pos_som_ent,
                             elt_val_som_ent,
                             NULL,
                             elt_val_color_ent,
                             cpt_coul_ent,
                             val_coul_ent,
                             cpt_elt_coul_ent);
}

/*----------------------------------------------------------------------------
 *  Lecture des groupes
 *----------------------------------------------------------------------------*/

static void
_ecs_pre_ideas__lit_groups(ecs_maillage_t   *maillage,
                           ecs_file_t       *fic_maillage,
                           int              *num_ligne,
                           unsigned          nbr_ent_par_ligne,
                           ecs_tab_int_t     tab_label_ent[],
                           bool              bool_label_ent_a_trier[])
{
  ecs_int_t    retour;                 /* Retour fonctions `scan' pour test   */
  char       chaine[ECS_LOC_LNG_MAX_CHAINE_IDEAS];             /* Ligne lue   */
  char     * ssch;                     /* Sous-chaine de lecture              */

  /* Variables Ideas lues */
  int          ideas_grp_num;          /* Numero Ideas du groupe              */
  unsigned     ideas_nbr_entities;     /* Nombre d'entites a lire par groupe  */
  char         ideas_grp_name[ECS_IDEAS_LEN_GROUP_NAME]; /* Nom du groupe     */
  int        * ideas_typ_entity;       /* Types Ideas des entites du groupe   */
  ecs_int_t  * ideas_tag_entity;       /* Etiquettes Ideas des entites du grp */

  ecs_descr_t  * descr_grp;            /* Pointeur sur descripteur de table   */
  unsigned     nbr_entities_lues;      /* Nbr d'entites du groupe lues        */
  unsigned     nbr_entities_ligne;     /* Nbr d'entites  par ligne lue        */
  ecs_int_t    num_ent;                /* Numero de l'entite                  */
  size_t       num_label;              /* Numero d'indice de l'etiquette      */
  long         tag_tmp;
  ecs_tab_int_t  tab_num   ;
  ecs_tab_int_t  tab_tag   ;
  ecs_tab_int_t  tab_label_elt;        /* Vecteur des labels des elements     */
                                       /*  pour l'ensemble des entites        */
  ecs_tab_int_t  tab_label_elt_trie ;
  ecs_tab_int_t  vect_renum_val_elt ;

  size_t       cpt_elt;
  size_t       nbr_elt;
  size_t       nbr_elt_ent ;
  size_t       ielt;                   /* Indice de boucle sur les elements   */
  int          ient;                   /* Indice de boucle sur les entites    */
  ecs_int_t    ient_sup;
  unsigned     ientity;                /* Indice de boucle sur les elements   */
                                       /*  d'un groupe IDEAS                  */
  size_t       ipos;
  size_t       ival;

  bool         bool_aff_grp;
  bool         bool_label_elt_a_trier ;
  size_t       ent_nbr_elt[ECS_N_ENTMAIL];
  ecs_int_t    min_val_ent[ECS_N_ENTMAIL];
  ecs_int_t    max_val_ent[ECS_N_ENTMAIL];
  ecs_int_t    tab_ent[ECS_N_ENTMAIL];
  ecs_int_t    ind_ent;
  ecs_int_t    min_ent = ECS_N_ENTMAIL;
  ecs_int_t    tmp_ent;

  size_t       ent_cpt_elt[ECS_N_ENTMAIL]; /* Nombre d'elements par entite  */

  /* Stockage avant transfert */
  ecs_int_t   *ent_val_grp[ECS_N_ENTMAIL]; /* Reference /entite et /groupe,
                                              les elts. appart. au groupe  */

  ecs_table_t  *table_grp = NULL;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Initialisations */
  /*=================*/

  for (ient = 0; ient < ECS_N_ENTMAIL; ient++)
    ent_nbr_elt[ient] = ecs_table__ret_elt_nbr(maillage->table_def[ient]);

  tab_tag.nbr = 1;
  tab_num.nbr = 1;
  ECS_MALLOC(tab_num.val, 1, ecs_int_t);

  /*-------------------------------------------------------------------------*/
  /* On cherche s'il est possible d'ordonner les labels des elements         */
  /*  sans effectuer de tri :                                                */
  /*  simplement en ordonnant les entites                                    */
  /*                                                                         */
  /* Ce n'est possible que si :                                              */
  /* - les labels pour chaque entite sont ordonnes                           */
  /*   (`bool_label_ent_a_trier[ient] == true')                          */
  /* - les labels entre les differentes entites ne se chevauchent pas        */
  /*                                                                         */
  /* Exemple :                                                               */
  /*  -    si les labels des aretes   sont compris entre 7 et  9             */
  /*  - et si les labels des faces    sont compris entre 3 et  5             */
  /*  - et si les labels des cellules sont compris entre 1 et  2             */
  /* en concatenant les labels des cellules, puis des faces, puis des aretes */
  /*  l'ensemble des labels sur les elements sera ordonne                    */
  /*-------------------------------------------------------------------------*/

  bool_label_elt_a_trier = false;

  nbr_elt = 0;
  for (ient = ECS_ENTMAIL_FAC; ient < ECS_N_ENTMAIL; ient++) {
      nbr_elt += tab_label_ent[ient].nbr;
      tab_ent[ient] = ient;
      if (bool_label_ent_a_trier[ient] == true)
        bool_label_elt_a_trier = true;
  }

  if (bool_label_elt_a_trier == false) {

    for (ient = ECS_ENTMAIL_CEL; ient >= ECS_ENTMAIL_FAC; ient--) {

      if (tab_label_ent[ient].nbr != 0) {

        min_val_ent[ient]= tab_label_ent[ient].val[0];
        max_val_ent[ient]= tab_label_ent[ient].val[tab_label_ent[ient].nbr - 1];

        min_ent = ient;

      }
    }

    for (ient = min_ent;
         ient < ECS_N_ENTMAIL && bool_label_elt_a_trier == false;
         ient++) {

      if (tab_label_ent[ient].nbr != 0) {

        for (ient_sup = ient + 1;
             ient_sup < ECS_N_ENTMAIL && bool_label_elt_a_trier == false;
             ient_sup++) {

          if (tab_label_ent[ient_sup].nbr != 0) {

            /* Si les valeurs des labels entre entites se chevauchent,
               il faut ordonner l'ensemble des labels sur les elements */

            if(! (   min_val_ent[tab_ent[ient_sup]] > max_val_ent[tab_ent[ient]]
                  || max_val_ent[tab_ent[ient_sup]] < min_val_ent[tab_ent[ient]])) {

              bool_label_elt_a_trier = true;

            }
            else if(max_val_ent[tab_ent[ient_sup]]<min_val_ent[tab_ent[ient]]){

              /* On echange l'ordre des entites    */
              /*  pour la concatenation des labels */

              tmp_ent           = tab_ent[ient];
              tab_ent[ient]     = tab_ent[ient_sup];
              tab_ent[ient_sup] = tmp_ent;

            }
            /* else : rien a faire (les entites sont bien ordonnees */

          } /* Fin : s'il y a des labels pour cette entite superieure */

        } /* Fin : boucle sur les entites superieures a l'entite courante */

      } /* Fin : s'il y a des labels pour cette entite */

    } /* Fin : boucle sur les entites */

  } /* Fin : si toute les entites ont des labels ordonnes */

  if (bool_label_elt_a_trier == true)
    for (ient = ECS_ENTMAIL_FAC; ient < ECS_N_ENTMAIL; ient++)
      tab_ent[ient] = ient;

  ECS_MALLOC(tab_label_elt.val, nbr_elt, ecs_int_t);
  tab_label_elt.nbr = nbr_elt;

  cpt_elt = 0;

  for (ient = ECS_ENTMAIL_FAC; ient < ECS_N_ENTMAIL; ient++) {

    for (ielt = 0; ielt <  tab_label_ent[tab_ent[ient]].nbr; ielt++)
      tab_label_elt.val[cpt_elt++] = tab_label_ent[tab_ent[ient]].val[ielt];

  }

  /*------------------------------------------------*/
  /* Ordination des elements suivant leur etiquette */
  /*------------------------------------------------*/

  ECS_MALLOC(vect_renum_val_elt.val, nbr_elt, ecs_int_t);
  vect_renum_val_elt.nbr = nbr_elt;

  if (bool_label_elt_a_trier == true) {

    tab_label_elt_trie = ecs_tab_int__trie_et_renvoie(tab_label_elt,
                                                      vect_renum_val_elt);
  }
  else {
    tab_label_elt_trie = tab_label_elt;
    for (ipos = 0; ipos < nbr_elt; ipos++)
      vect_renum_val_elt.val[ipos] = ipos;
  }

  ideas_typ_entity   = NULL;  /* Pour le REALLOC au 1er passage */
  ideas_tag_entity   = NULL;  /* Pour le REALLOC au 1er passage */
  ideas_nbr_entities = 0;

  /*========================================================*/
  /* Tant qu'on n'a pas atteint la fin de fichier           */
  /* et que la ligne lue n'est pas le separateur de dataset */
  /*========================================================*/

  while (ecs_file_gets_try(chaine, ECS_LOC_LNG_MAX_CHAINE_IDEAS,
                           fic_maillage, num_ligne) != NULL
         && strcmp(chaine, ECS_IDEAS_SEPARATEUR_DATASET) != 0) {

    /* On alloue et initialise pour le groupe a lire */
    /*===============================================*/

    for (ient = 0; ient < ECS_N_ENTMAIL; ient++) {

      ent_cpt_elt[ient] = 0;
      ent_val_grp[ient] = NULL;

      if (ent_nbr_elt[ient] != 0)
        ECS_MALLOC(ent_val_grp[ient], ent_nbr_elt[ient], ecs_int_t);

      for (ival = 0; ival < ent_nbr_elt[ient]; ival++)
        ent_val_grp[ient][ival] = 0;
    }

    /* lecture d'un groupe */
    /*=====================*/

    /* Lecture des caracteristiques du groupe */
    /*----------------------------------------*/

    /* Format Ideas : 8I10 */

    retour = sscanf(chaine,
                    " %d %*d %*d %*d %*d %*d %*d %u",
                    &ideas_grp_num,
                    &ideas_nbr_entities);

    if (retour != 2)
        ecs_error(__FILE__, __LINE__, errno,
                  _("Error reading line %d of file \"%s\"."),
                  *num_ligne, ecs_file_get_name(fic_maillage));

    /* Lecture du nom du groupe */
    /*--------------------------*/

    /* Format Ideas : 20A2 */

    ecs_file_gets(chaine, ECS_LOC_LNG_MAX_CHAINE_IDEAS,
                  fic_maillage, num_ligne);

    /* retour = sscanf(chaine, "%[^\n\f\r]",ideas_grp_name); */

    strcpy(ideas_grp_name, "\0");
    ssch = strtok(chaine, " \n");
    strcat(ideas_grp_name, ssch);
    while ((ssch = strtok(NULL, " \n")) != NULL) {
      strcat(ideas_grp_name, " ");
      strcat(ideas_grp_name, ssch);
    }

    /* Lecture des caracteristiques des elements du groupe */
    /*-----------------------------------------------------*/

    /* Format Ideas : 8I10 (pour chaque ligne) */

    ECS_REALLOC(ideas_typ_entity, ideas_nbr_entities, int);
    ECS_REALLOC(ideas_tag_entity, ideas_nbr_entities, ecs_int_t);

    nbr_entities_lues = 0;

    while (nbr_entities_lues != ideas_nbr_entities) {

      ecs_file_gets(chaine, ECS_LOC_LNG_MAX_CHAINE_IDEAS,
                    fic_maillage, num_ligne);

      /* Lecture des types et numeros des elements pour la ligne courante */

      nbr_entities_ligne
        = ECS_MIN(ideas_nbr_entities - nbr_entities_lues, nbr_ent_par_ligne);

      ssch = strtok(chaine, " ");

      for (ientity = 0; ientity < nbr_entities_ligne; ientity++) {

        retour = sscanf(ssch, "%d",
                        &ideas_typ_entity[nbr_entities_lues + ientity]);
        ssch   = strtok(NULL, " ");
        retour = sscanf(ssch, "%ld", &tag_tmp);
        ideas_tag_entity[nbr_entities_lues + ientity] = tag_tmp;
        ssch   = strtok(NULL, " ");

        if (nbr_ent_par_ligne == ECS_IDEAS_NBR_GRP_ENT_PER_LINE2) {

          retour = sscanf(ssch, "%*d");
          ssch   = strtok(NULL, " ");
          retour = sscanf(ssch, "%*d");
          ssch   = strtok(NULL, " ");

        }
      }

      if (ssch != NULL)
        ecs_error(__FILE__, __LINE__, errno,
                  _("Error reading line %d of file \"%s\"."),
                  *num_ligne, ecs_file_get_name(fic_maillage));

      nbr_entities_lues += nbr_entities_ligne;

    }

    /* Determination de l'entite a laquelle correspondent */
    /*  les elements du groupe                            */
    /*----------------------------------------------------*/

    for (ientity = 0; ientity < ideas_nbr_entities; ientity++) {

      tab_tag.val = &ideas_tag_entity[ientity];

      if (ideas_typ_entity[ientity] == ECS_IDEAS_TYP_CODE_ELEMENTS) {

        /*---------------------------------------*/
        /* Recherche de l'indice correspondant   */
        /*  au label `ideas_typ_entity[ientity]' */
        /*---------------------------------------*/

        ecs_tab_int__recherche(tab_tag,
                               tab_label_elt_trie,
                               tab_num);

        if (tab_num.val[0] == -1)
          continue;

        num_label = vect_renum_val_elt.val[tab_num.val[0]];

        ind_ent = ECS_ENTMAIL_FAC;
        num_ent = tab_ent[ind_ent];
        nbr_elt_ent = ent_nbr_elt[num_ent];
        while (   ind_ent   < ECS_N_ENTMAIL
               && (num_label + 1) > nbr_elt_ent) {
          num_label -=  nbr_elt_ent;
          num_ent = tab_ent[++ind_ent];
          nbr_elt_ent = ent_nbr_elt[num_ent];
        }

        assert(ind_ent != ECS_N_ENTMAIL);

        /* Stockage des valeurs lues avant transfert dans maillage */

        ent_val_grp[num_ent][num_label] = 1;

        /* Incrementation du nombre d'objets lus */

        ent_cpt_elt[num_ent]++;

      }

    } /* Fin : boucle sur les entites Ideas du groupe */

    /* Boucle de remplissage des entites du maillage */
    /*===============================================*/

    bool_aff_grp = false;

    for (ient = ECS_ENTMAIL_FAC; ient < ECS_N_ENTMAIL; ient++) {

      /*--------------------------------------------------------------------*/
      /* S'il y a au moins un element du groupe de cette entite de maillage */
      /*--------------------------------------------------------------------*/

      if (ent_cpt_elt[ient] != 0) {

        bool_aff_grp = true;

        assert(ent_cpt_elt[ient] <= ent_nbr_elt[ient]);

        /* Creation du descripteur de table correspondant au groupe lu */
        /*-------------------------------------------------------------*/

        descr_grp = ecs_descr__cree(ECS_DESCR_IDE_NUL,
                                    ideas_grp_name);

        /* Transformation du tableau referencant le groupe en une table */
        /*--------------------------------------------------------------*/

        table_grp = ecs_table__transforme_tableau(ent_nbr_elt[ient],
                                                  ent_val_grp[ient],
                                                  descr_grp);

        if (maillage->table_att[ient] == NULL)
          maillage->table_att[ient] = table_grp;
        else
          ecs_table_att__assemble(maillage->table_att[ient],
                                  table_grp);

      } /* Fin si le nombre d'elements referencant le groupe n'est pas nul */

    } /* Fin de la boucle sur les entites */

    /* Affichage du bilan des donnees lues pour les groupes */
    /*======================================================*/

    if (bool_aff_grp == true)
      printf("  %s %d \"%s\"\n",
             _("Group"), ideas_grp_num, ideas_grp_name);

    ecs_maillage_pre__aff_nbr_par_ent(0,
                                      ent_cpt_elt,
                                      0);

    /* Incrementation du compteur sur les groupes */
    /*--------------------------------------------*/

    for (ient = 0; ient < ECS_N_ENTMAIL; ient++)
      if (ent_nbr_elt[ient] != 0)
        ECS_FREE(ent_val_grp[ient]);

    /* Re-initialisation des compteurs par entite pour le groupe suivant */
    for (ient = 0; ient < ECS_N_ENTMAIL; ient++)
      ent_cpt_elt[ient] = 0;


  } /* Fin : tant qu'il y a une ligne a lire et                       */
    /*       tant que la ligne lue n'est pas le separateur de dataset */

  if (strcmp(chaine, ECS_IDEAS_SEPARATEUR_DATASET) != 0)
    ecs_file_read_check_error(fic_maillage, *num_ligne);

  /* Liberation des tableaux locaux */
  /*================================*/

  ECS_FREE(tab_label_elt.val);
  ECS_FREE(tab_num.val);

  if (bool_label_elt_a_trier == true)
    ECS_FREE(tab_label_elt_trie.val);

  ECS_FREE(vect_renum_val_elt.val);

  ECS_FREE(ideas_typ_entity);
  ECS_FREE(ideas_tag_entity);
}

/*============================================================================
 *  Fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Lecture d'un fichier I-DEAS Master Series au format universel
 *   et affectation des donnees dans la structure de maillage
 *
 *  Hypothèses de lecture :
 *   le dataset sur les groupes doit se trouver placer dans le fichier Ideas
 *   après les dataset sur les noeuds et sur les éléments ; le dataset
 *   sur noeuds doit lui-même se trouver après celui sur les systèmes de
 *   coordonnées.
 *----------------------------------------------------------------------------*/

ecs_maillage_t *
ecs_pre_ideas__lit_maillage(const char  *nom_fic_maillage)
{
  ecs_file_t  *fic_maillage;
  ecs_int_t    retour;                   /* Retour fnctions `scan' pour test  */
  char         chaine[ECS_LOC_LNG_MAX_CHAINE_IDEAS];           /* Ligne lue   */
  bool         bool_dataset_suivant;     /* Indicateur de lecture de dataset  */
  bool         bool_dataset_elements;    /* Indicateur de lecture du dataset
                                            sur les elements */
  bool         bool_dataset_noeuds;      /* Indicateur de lecture du dataset
                                            sur les noeuds */
  ecs_int_t    ient;                     /* Indice de boucle sur les entites */
  unsigned     entities_per_line;        /* Nbr d'entites Ideas par ligne */
  int          num_ligne;                /* Compteur des lignes lues */
  int          num_dataset;              /* Numero du dataset lu */
  bool         bool_label_som_a_trier;
  bool         bool_label_ent_a_trier[ECS_N_ENTMAIL];     /* Ind labels tries */
  ecs_tab_int_t  tab_label_ent[ECS_N_ENTMAIL];    /* Tabl. labels par entite */
  ecs_tab_int_t  tab_ind_coo_sys;        /* Corresp. label -> sys. coo référ. */
  ecs_int_t      *som_val_label;
  _ecs_ideas_sys_coord_t  *coo_sys_ref;  /* Syst. coord. de référence */

  /* Création d'un maillage initialement vide (valeur de retour) */

  ecs_maillage_t  *maillage = ecs_maillage__cree_nodal();

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Affichage du titre */
  /*====================*/

  printf(_("\n\n"
           "Reading mesh from file in I-deas universal format\n"
           "----------------------\n"));

  printf(_("  Mesh file: %s\n\n\n"), nom_fic_maillage);

  /* Initialisations */
  /*=================*/

  bool_label_som_a_trier = false;

  for (ient = 0; ient < ECS_N_ENTMAIL; ient++) {
    bool_label_ent_a_trier[ient] = false;
    tab_label_ent[ient].val = NULL;
    tab_label_ent[ient].nbr = 0;
  }

  tab_ind_coo_sys.nbr = 0;
  tab_ind_coo_sys.val = NULL;
  coo_sys_ref = NULL;

  num_ligne = 1;

  som_val_label = NULL;

  /* Ouverture du fichier Ideas en lecture */
  /*---------------------------------------*/

  fic_maillage = ecs_file_open(nom_fic_maillage,
                               ECS_FILE_MODE_READ,
                               ECS_FILE_TYPE_TEXT);

  bool_dataset_suivant  = false;

  bool_dataset_elements = false;
  bool_dataset_noeuds   = false;

  /*================================================*/
  /* Boucle sur les lignes du fichier de maillage : */
  /* tant qu'on n'a pas atteint la fin de fichier   */
  /*================================================*/

  while (ecs_file_gets_try(chaine, ECS_LOC_LNG_MAX_CHAINE_IDEAS,
                           fic_maillage, &num_ligne) != NULL) {

    /* Si la chaine lue est la chaine de separation des dataset de Ideas */
    /*===================================================================*/

    if ( strcmp(chaine, ECS_IDEAS_SEPARATEUR_DATASET) == 0) {

      if (bool_dataset_suivant == true) {
        /* On a lu le separateur de fin du dataset precedent           */
        /* On lit maintenant le separateur de debut du dataset suivant */
        if (ecs_file_gets_try(chaine, ECS_LOC_LNG_MAX_CHAINE_IDEAS,
                              fic_maillage, &num_ligne) == NULL)
          bool_dataset_suivant = false;
      }
      else {
        /* c'est le premier separateur lu : celui du debut de fichier */
        bool_dataset_suivant = true;
      }

      if (bool_dataset_suivant == true) {

        /* Lecture du numero du dataset */
        /*------------------------------*/

        ecs_file_gets(chaine, ECS_LOC_LNG_MAX_CHAINE_IDEAS,
                      fic_maillage, &num_ligne);

        retour = sscanf(chaine, " %d", &num_dataset);

        if ( retour != 1)
          ecs_error(__FILE__, __LINE__, errno,
                    _("Error reading line %d of file \"%s\"."),
                    num_ligne, ecs_file_get_name(fic_maillage));

        /* Suivant le numero du dataset lu ... */
        /*-------------------------------------*/

        switch(num_dataset) {

        case ECS_IDEAS_DATASET_SYS_COORD_2420:

          /* Lecture de la définition des systèmes de coordonnées */

          _ecs_pre_ideas__lit_sys_coord(fic_maillage,
                                        &num_ligne,
                                        &tab_ind_coo_sys,
                                        &coo_sys_ref);

          /* On a déjà lu dans la fonction le separateur de fin de dataset */
          bool_dataset_suivant = false;

          break;

        case ECS_IDEAS_DATASET_NODES_2411:

          /* Lecture des coordonnees des noeuds */
          /*------------------------------------*/

          _ecs_pre_ideas__lit_nodes(maillage,
                                    fic_maillage,
                                    &num_ligne,
                                    &som_val_label,
                                    &bool_label_som_a_trier,
                                    tab_ind_coo_sys,
                                    coo_sys_ref);

          /* On a deja lu dans la fonction le separateur de fin de dataset */
          bool_dataset_suivant = false;

          /* On indique qu'on a lu le dataset sur les noeuds */
          bool_dataset_noeuds  = true;

          break;

        case ECS_IDEAS_DATASET_ELEMENTS_2412:

          /* Lecture de la table de connectivite */
          /*-------------------------------------*/

          _ecs_pre_ideas__lit_elements(maillage,
                                       fic_maillage,
                                       &num_ligne,
                                       maillage->n_vertices,
                                       &som_val_label,
                                       tab_label_ent,
                                       bool_label_ent_a_trier);

          /* On a deja lu dans la fonction le separateur de fin de dataset */
          bool_dataset_suivant = false;

          /* On indique qu'on a lu le dataset sur les elements */
          bool_dataset_elements = true;

          break;

        case ECS_IDEAS_DATASET_GROUPS_2430:
        case ECS_IDEAS_DATASET_GROUPS_2432:
        case ECS_IDEAS_DATASET_GROUPS_2435:
        case ECS_IDEAS_DATASET_GROUPS_2452:
        case ECS_IDEAS_DATASET_GROUPS_2467:
        case ECS_IDEAS_DATASET_GROUPS_2477:

          /* La construction du maillage fait l'hypothese                */
          /*  que les noeuds et les elements sont lus avant les groupes. */
          /* On verifie qu'on a lu les datasets sur les noeuds et        */
          /* les elements avant de lire le dataset sur les groupes       */
          /*-------------------------------------------------------------*/

          if (bool_dataset_noeuds == false ||
              bool_dataset_elements == false) {

            ecs_error(__FILE__, __LINE__, 0,
                      _("Error reading an I-deas universal file:\n"
                        "at line %d of file \"%s\".\n"
                        "Node and element datasets must be read\n"
                        "before the groups dataset."),
                      (int)num_ligne, ecs_file_get_name(fic_maillage));
          }

          /* Lecture des groupes dans un format obsolete */
          /*                  ou dans le format actuel   */
          /*---------------------------------------------*/

          if (   num_dataset == ECS_IDEAS_DATASET_GROUPS_2430
              || num_dataset == ECS_IDEAS_DATASET_GROUPS_2432) {
            /* Lecture des groupes dans un format obsolete */
            entities_per_line = ECS_IDEAS_NBR_GRP_ENT_PER_LINE4;
          }
          else {
            /* Lecture des groupes dans le format actuel   */
            entities_per_line = ECS_IDEAS_NBR_GRP_ENT_PER_LINE2;
          }

          _ecs_pre_ideas__lit_groups(maillage,
                                     fic_maillage,
                                     &num_ligne,
                                     entities_per_line,
                                     tab_label_ent,
                                     bool_label_ent_a_trier);

          /* On a deja lu dans la fonction le separateur de fin de dataset */
          bool_dataset_suivant = false;

          break;

        default:

          ;  /* On poursuit la lecture */

        } /* Fin : `switch(num_dataset)' : numero du dataset lu ... */


      } /* Fin : si on a lu le separateur de dataset */

    } /* Fin : on n'a pas atteint la fin de fichier au cours de la boucle */

  } /* Fin : boucle tant qu'on n'a pas atteint la fin de fichier  */

  if (ecs_file_eof(fic_maillage) == 0)
    ecs_error(__FILE__, __LINE__, errno,
              _("Error reading line %d of file \"%s\"."),
              num_ligne, ecs_file_get_name(fic_maillage));

  /* else : la fin de fichier a bien ete atteinte */

  /* On verifie qu'on a bien lu des noeuds et des elements */

  if (bool_dataset_noeuds == false)
    ecs_error(__FILE__, __LINE__, 0,
              _("Error reading an I-deas universal file:\n"
                "at line %d of file \"%s\".\n"
                "Dataset \"%d\" containing the nodes "
                "definition has not been found."),
              (int)num_ligne, ecs_file_get_name(fic_maillage),
              ECS_IDEAS_DATASET_NODES_2411);

  if (bool_dataset_elements == false)
    ecs_error(__FILE__, __LINE__, 0,
              _("Error reading an I-deas universal file:\n"
                "at line %d of file \"%s\".\n"
                "Dataset \"%d\" containing the elements "
                "definition has not been found."),
              num_ligne, ecs_file_get_name(fic_maillage),
              ECS_IDEAS_DATASET_ELEMENTS_2412);

  /* Fermeture du fichier de lecture du maillage */
  /*---------------------------------------------*/

  ecs_file_free(fic_maillage);

  /* Liberations */
  /*=============*/

  tab_ind_coo_sys.nbr = 0;
  ECS_FREE(tab_ind_coo_sys.val);
  ECS_FREE(coo_sys_ref);

  for (ient = 0; ient < ECS_N_ENTMAIL; ient++)
    if (tab_label_ent[ient].val != NULL)
      ECS_FREE(tab_label_ent[ient].val);

  return maillage;
}

/*----------------------------------------------------------------------------*/
