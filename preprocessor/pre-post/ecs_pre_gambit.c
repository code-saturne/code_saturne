/*============================================================================
 *  Définition de la fonction
 *   de lecture d'un fichier de maillage GAMBIT neutral
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

#include "ecs_pre_gambit.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' privés   du  paquetage courant
 *----------------------------------------------------------------------------*/


/*============================================================================
 *  Structures locales
 *============================================================================*/

/*============================================================================
 *                       Définition de constantes globales statiques
 *============================================================================*/


/* Liste des identificateurs de description des elements finis */
/*=============================================================*/

#define ECS_GAMBIT_EDGE            1
#define ECS_GAMBIT_QUADRILATERAL   2
#define ECS_GAMBIT_TRIANGLE        3
#define ECS_GAMBIT_BRICK           4
#define ECS_GAMBIT_WEDGE           5
#define ECS_GAMBIT_TETRAHEDRON     6
#define ECS_GAMBIT_PYRAMID         7

/* Définition des éléments */
/*=========================*/

typedef struct {

  ecs_int_t       gambit_typ   ; /* Type GAMBIT de l'élément  */
  ecs_elt_typ_t   ecs_typ      ; /* Type ECS    de l'élément  */
  ecs_int_t       nbr_som[4]   ; /* Nombre de sommets GAMBIT (variantes 1-4) */
  ecs_int_t       ind_som[4][8]; /* Numéros de sommets ECS   (variantes 1-4) */
  ecs_int_t       nbr_sselt    ; /* Nombre de sous-éléments */
  ecs_sous_elt_t  sous_elt[6]  ; /* Déf. sous-éléments (numérotation locale
                                   * usuelle de l'Enveloppe, mais de 0 à n-1
                                   * au lieu de 1 à n) */

} ecs_loc_gambit_elt_t;


static const ecs_loc_gambit_elt_t  ecs_loc_gambit_elt_liste_c[7] = {

  {                        /* 1 */
    ECS_GAMBIT_EDGE,
    ECS_ELT_TYP_NUL,
    { 2, 3, 0, 0 },
    {
      { 0 },
      { 0 }
    },
    0,
    {
      {ECS_ELT_TYP_NUL  , {0}}
    }
  },
  {                        /* 2 */
    ECS_GAMBIT_QUADRILATERAL,
    ECS_ELT_TYP_FAC_QUAD,
    { 4, 8, 9, 0 },
    {
      { 0, 1, 2, 3 },
      { 0, 2, 4, 6 },
      { 0, 2, 4, 6 }                               /*  3 x-------x 2          */
    },                                             /*    |       |            */
    0,                                             /*    |       |            */
    {                                              /*    |       |            */
      {ECS_ELT_TYP_NUL  , {0}}                     /*  0 x-------x 1          */
    }
  },
  {                        /* 3 */
    ECS_GAMBIT_TRIANGLE,
    ECS_ELT_TYP_FAC_TRIA,
    { 3, 6, 7, 0 },
    {
      { 0, 1, 2 },
      { 0, 2, 4 },
      { 0, 2, 4 }                                  /*        x 2              */
    },                                             /*       / \               */
    0,                                             /*      /   \              */
    {                                              /*     /     \             */
      {ECS_ELT_TYP_NUL  , {0}}                     /*  0 x-------x 1          */
    }
  },
  {                        /* 4 */
    ECS_GAMBIT_BRICK,
    ECS_ELT_TYP_CEL_HEXA,
    { 8, 20, 27, 0 },
    {
      {  4,  5,  1,  0,  6,  7,  3,  2 },
      { 12, 14,  2,  0, 17, 19,  7,  5 },
      { 18, 20,  2,  0, 24, 26,  8,  6 }
    },
    6,                                             /*     7 x-------x 6       */
    {                                              /*      /|      /|         */
      {ECS_ELT_TYP_FAC_QUAD, { 3, 2, 1, 0 }},      /*     / |     / |         */
      {ECS_ELT_TYP_FAC_QUAD, { 2, 6, 5, 1 }},      /*  4 x-------x5 |         */
      {ECS_ELT_TYP_FAC_QUAD, { 6, 7, 4, 5 }},      /*    | 3x----|--x 2       */
      {ECS_ELT_TYP_FAC_QUAD, { 7, 3, 0, 4 }},      /*    | /     | /          */
      {ECS_ELT_TYP_FAC_QUAD, { 2, 3, 7, 6 }},      /*    |/      |/           */
      {ECS_ELT_TYP_FAC_QUAD, { 0, 1, 5, 4 }}       /*  0 x-------x 1          */
    }
  },
  {                       /*  5 */
    ECS_GAMBIT_WEDGE,
    ECS_ELT_TYP_CEL_PRISM,
    { 6, 15, 18, 0 },
    {
      { 0,  1,  2,  3,  4,  5 },
      { 0,  2,  5,  9, 11, 14 },
      { 0,  2,  5, 12, 14, 17 }
    },                                             /*  3 x-------x 5          */
    5,                                             /*    |\     /|            */
    {                                              /*    | \   / |            */
      {ECS_ELT_TYP_FAC_QUAD, { 0, 1, 4, 3 }},      /*  0 x- \-/ -x 2          */
      {ECS_ELT_TYP_FAC_QUAD, { 1, 2, 5, 4 }},      /*     \ 4x  /             */
      {ECS_ELT_TYP_FAC_QUAD, { 2, 0, 3, 5 }},      /*      \ | /              */
      {ECS_ELT_TYP_FAC_TRIA, { 0, 2, 1 }   },      /*       \|/               */
      {ECS_ELT_TYP_FAC_TRIA, { 3, 4, 5 }   }       /*        x 1              */
    }
  },
  {                        /* 6 */
    ECS_GAMBIT_TETRAHEDRON,
    ECS_ELT_TYP_CEL_TETRA,
    { 4, 10, 0, 0 },
    {
      { 0,  1,  2,  3 },
      { 0,  2,  5,  9 }
    },                                             /*        x 3              */
    4,                                             /*       /|\               */
    {                                              /*      / | \              */
      {ECS_ELT_TYP_FAC_TRIA, { 1, 0, 2 }},         /*     /  |  \             */
      {ECS_ELT_TYP_FAC_TRIA, { 0, 1, 3 }},         /*  0 x- -|- -x 2          */
      {ECS_ELT_TYP_FAC_TRIA, { 1, 2, 3 }},         /*     \  |  /             */
      {ECS_ELT_TYP_FAC_TRIA, { 2, 0, 3 }}          /*      \ | /              */
    }                                              /*       \|/               */
  },                                               /*        x 1              */
  {                        /* 7 */
    ECS_GAMBIT_PYRAMID,
    ECS_ELT_TYP_CEL_PYRAM,
    { 5, 13, 14, 18 },
    {
      { 0,  1,  3,  2,  4 },
      { 0,  2,  7,  5, 12 },
      { 0,  2,  8,  6, 12 },
      { 0,  2,  8,  6, 17 }
    },                                             /*         4 x             */
    5,                                             /*          /|\            */
    {                                              /*         //| \           */
      {ECS_ELT_TYP_FAC_QUAD, { 0, 3, 2, 1 }},      /*        // |  \          */
      {ECS_ELT_TYP_FAC_TRIA, { 0, 1, 4 }   },      /*     3 x/--|---x 2       */
      {ECS_ELT_TYP_FAC_TRIA, { 1, 2, 4 }   },      /*      //   |  /          */
      {ECS_ELT_TYP_FAC_TRIA, { 2, 3, 4 }   },      /*     //    | /           */
      {ECS_ELT_TYP_FAC_TRIA, { 3, 0, 4 }   }       /*  0 x-------x 1          */
    },
  }
};


/* Definition des groupes */
/*========================*/

        /* Longueur des noms de groupe                    */
#define ECS_GAMBIT_LEN_GROUP_NAME                                    40

        /* Nombre max d'entites par ligne                 */
        /* pour la definition des groupes                 */
#define ECS_GAMBIT_NBR_GRP_ENT_PER_LINE2                              2
        /* Idem pour les anciens dataset                  */
#define ECS_GAMBIT_NBR_GRP_ENT_PER_LINE4                              4


/* Codes des entites appartenant a un groupe */
/*-------------------------------------------*/

#define ECS_GAMBIT_TYP_CODE_NODES                                     7
#define ECS_GAMBIT_TYP_CODE_ELEMENTS                                  8


/*============================================================================
 *  Définitions de parametres-macros
 *============================================================================*/

/* Pour une lecture de 80 caracteres par ligne  */
/* auxquels il faut ajouter le `\n' et le `\0'  */
/* pour l'affectation dans la chaine receptrice */
/* On ajoute 2 caracteres de securite pour les  */
/* fichiers au format 'DOS'                     */
#define ECS_LOC_LNG_MAX_CHAINE_GAMBIT  84        /* Dimension des chaines */

/*============================================================================
 *  Fonctions privées
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Remplissage d'une chaîne de taille ECS_LOC_LNG_MAX_CHAINE_GAMBIT
 *  par des blancs (initialisation)
 *----------------------------------------------------------------------------*/

static void ecs_loc_pre_gambit__raz_chaine
(
 char                         *const chaine
)
{
  int ind;

  for (ind = 0; ind < ECS_LOC_LNG_MAX_CHAINE_GAMBIT; ind++)
    chaine[ind] = ' ';

  chaine[ECS_LOC_LNG_MAX_CHAINE_GAMBIT - 1] = '\0';
}

/*----------------------------------------------------------------------------
 *  Lecture de la fin d'une rubrique
 *----------------------------------------------------------------------------*/

static void ecs_loc_pre_gambit__fin_section
(
 ecs_file_t                   *const fic_maillage , /* --> Descr. fichier     */
 int                          *const num_ligne      /* <-> Cpt. lignes lues   */
)
{
  char  chaine[ECS_LOC_LNG_MAX_CHAINE_GAMBIT];

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


  ecs_file_gets(chaine, ECS_LOC_LNG_MAX_CHAINE_GAMBIT,
                fic_maillage, num_ligne);

  if (strncmp(chaine, "ENDOFSECTION", 12) != 0)
    ecs_error(__FILE__, __LINE__, 0,
              _("Error reading line %d of file \"%s\";\n"
                "An end of section (ENDOFSECTION) was expected."),
              *num_ligne, ecs_file_get_name(fic_maillage));

}


/*----------------------------------------------------------------------------
 *  Saut d'une rubrique
 *----------------------------------------------------------------------------*/

static void ecs_loc_pre_gambit__saut_section
(
 ecs_file_t                   *const fic_maillage , /* --> Descr. fichier     */
 int                          *const num_ligne      /* <-> Cpt. lignes lues   */
)
{
  char  chaine[ECS_LOC_LNG_MAX_CHAINE_GAMBIT];

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  do {
    ecs_file_gets(chaine, ECS_LOC_LNG_MAX_CHAINE_GAMBIT,
                  fic_maillage, num_ligne);
  } while (strncmp(chaine, "ENDOFSECTION", 12) != 0);

}


/*----------------------------------------------------------------------------
 *  Lecture de l'entête
 *----------------------------------------------------------------------------*/

static void ecs_loc_pre_gambit__lit_entete
(
 ecs_file_t                   *const fic_maillage , /* --> Descr. fichier     */
 int                          *const num_ligne    , /* <-> Cpt. lignes lues   */
 int                          *const dim_e        , /* --> Dim. espace        */
 ecs_int_t                    *const nbr_som      , /* --> Nb. sommets        */
 ecs_int_t                    *const nbr_elt      , /* --> Nb. éléments       */
 ecs_int_t                    *const nbr_grp      , /* --> Nb. groupes        */
 ecs_int_t                    *const nbr_cl         /* --> Nb. C.L.           */
)
{
  ecs_int_t  retour;                 /* Retour fonctions `scan' pour test    */
  char       chaine[ECS_LOC_LNG_MAX_CHAINE_GAMBIT];           /* Ligne lue   */
  char       sch1[ECS_LOC_LNG_MAX_CHAINE_GAMBIT];
  char       sch2[ECS_LOC_LNG_MAX_CHAINE_GAMBIT];
  ecs_int_t  ind;

  /* Variables GAMBIT lues */

  int numnp, nelem, ngprs, nbsets, ndfcd, ndfvl;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


  /* Format */
  /*--------*/

  ecs_file_gets(chaine, ECS_LOC_LNG_MAX_CHAINE_GAMBIT,
                fic_maillage, num_ligne);

  if (strncmp(chaine, "** GAMBIT NEUTRAL FILE", 22) != 0)
    ecs_error(__FILE__, __LINE__, 0,
              _("Format error for file \"%s\":\n"
                "This file does not seem to be in GAMBIT neutral format\n"
                "(line 2 does not start with \"** GAMBIT NEUTRAL FILE\")."),
                ecs_file_get_name(fic_maillage));

  /* Titre */
  /*-------*/

  ecs_loc_pre_gambit__raz_chaine(chaine);

  ecs_file_gets(chaine, ECS_LOC_LNG_MAX_CHAINE_GAMBIT,
                fic_maillage, num_ligne);

  /*
    Suppression blancs en fin de titre (80 caractères en général rarement
    utilisés) pour affichage sur ligne plus courte
  */

  chaine[80] = '\0';
  for (ind = 80; ind > 0 && chaine[ind - 1] == ' '; ind--)
    chaine[ind] = '\0';

  printf(_("  Title        : %.80s\n"), chaine);


  /* Source */
  /*--------*/

  ecs_file_gets(chaine, ECS_LOC_LNG_MAX_CHAINE_GAMBIT,
                fic_maillage, num_ligne);

  retour = sscanf(chaine, "%*s %s %*s %s", sch1, sch2);

  if (retour == 2)
    printf(_("  Created with : %s %s\n"), sch1, sch2);


  /* Date et heure */
  /*---------------*/

  ecs_file_gets(chaine, ECS_LOC_LNG_MAX_CHAINE_GAMBIT,
                fic_maillage, num_ligne);

  chaine[80] = '\0';
  for (ind = 80; ind > 0 && chaine[ind - 1] == ' '; ind--)
    chaine[ind] = '\0';
  for (ind = 0; chaine[ind] != '\0' && chaine[ind] == ' '; ind++);

  printf(_("  Date         : %s\n"), chaine + ind);


  /* Infos sur la taille */
  /*---------------------*/

  /* Ligne "     NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL" */

  ecs_file_gets(chaine, ECS_LOC_LNG_MAX_CHAINE_GAMBIT,
                fic_maillage, num_ligne);

  /* Ligne avec valeurs à lire */

  ecs_file_gets(chaine, ECS_LOC_LNG_MAX_CHAINE_GAMBIT,
                fic_maillage, num_ligne);

  retour = sscanf(chaine, "%d %d %d %d %d %d",
                  &numnp, &nelem, &ngprs, &nbsets, &ndfcd, &ndfvl);

  if (retour != 6)
    ecs_error(__FILE__, __LINE__, 0,
              _("Error reading line %d of file \"%s\"."),
              *num_ligne, ecs_file_get_name(fic_maillage));

  printf(_("  Initial data: %10d points\n"
           "                %10d elements\n"
           "                %10d groups\n"
           "                %10d boundary conditions\n"
           "                spatial dimension: %d\n"
           "                velocity components: %d\n\n"),
         numnp, nelem, ngprs, nbsets, ndfcd, ndfvl);

  /* Fin */
  /*-----*/

  ecs_loc_pre_gambit__fin_section(fic_maillage, num_ligne);

  *dim_e   = ndfcd;
  *nbr_som = numnp;
  *nbr_elt = nelem;
  *nbr_grp = ngprs;
  *nbr_cl  = nbsets;

}


/*----------------------------------------------------------------------------
 *  Lecture des coordonnées des sommets
 *----------------------------------------------------------------------------*/

static void
ecs_loc_pre_gambit__lit_coords(ecs_file_t    *fic_maillage,
                               int            dim_e,
                               ecs_int_t      nbr_som,
                               ecs_coord_t    som_val_coord[],
                               ecs_int_t      som_val_label[],
                               int           *num_ligne)
{
  ecs_int_t    isom;
  ecs_int_t    icoo;

  /* Variables Gambit lues */

  char         chaine[ECS_LOC_LNG_MAX_CHAINE_GAMBIT];
  int          num_nod;
  double       coord[3];


  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


  /*=========================================================*/
  /* Boucle de lecture des numéros et coordonnées des points */
  /*=========================================================*/

  for (isom = 0; isom < nbr_som; isom++) {

    ecs_int_t  retour = 0;

    ecs_file_gets(chaine, ECS_LOC_LNG_MAX_CHAINE_GAMBIT,
                  fic_maillage, num_ligne);

    if (dim_e == 3)
      retour = sscanf(chaine, "%d %le %le %le",
                      &num_nod, coord, coord+1, coord+2);

    else if (dim_e == 2) {
      retour = sscanf(chaine, "%d %le %le",
                      &num_nod, coord, coord+1);
      coord[2] = 0.0;
    }

    if (retour != (int)(dim_e + 1))
      ecs_error(__FILE__, __LINE__, 0,
                _("Error reading line %d of file \"%s\"\n"
                  "while decoding point %d (of %d)."),
                *num_ligne, ecs_file_get_name(fic_maillage),
                (int)(isom + 1), (int)nbr_som);

    /* Coordonnees du sommet lu */

    for (icoo = 0; icoo < 3; icoo++)
      som_val_coord[isom * 3 + icoo] = coord[icoo];


    /* Etiquette   du sommet lu */

    som_val_label[isom] = num_nod;

  }

  /* Fin de la rubrique */

  ecs_loc_pre_gambit__fin_section(fic_maillage, num_ligne);

  /* Maillage en dimension 3 */

  if (dim_e == 2) {
    ecs_warn();
    printf(_("The mesh is 2d."));
  }

}


/*----------------------------------------------------------------------------
 *  Lecture de la table de connectivité
 *----------------------------------------------------------------------------*/

static void
ecs_loc_pre_gambit__lit_elements(ecs_file_t  *fic_maillage,
                                 size_t       nbr_elt,
                                 int         *num_ligne,
                                 size_t       cpt_elt_ent[],
                                 ecs_size_t  *elt_pos_som_ent[],
                                 ecs_int_t   *elt_val_som_ent[],
                                 ecs_int_t   *elt_val_typ_geo_ent[],
                                 ecs_int_t   *elt_val_label_ent[])
{
  ecs_int_t    retour;               /* Retour fonctions `scan' pour test    */
  char         chaine[ECS_LOC_LNG_MAX_CHAINE_GAMBIT];          /* Ligne lue  */

  /* Variables GAMBIT lues */
  int          ne;                   /* Numéro (label) de l'élément          */
  int          ntype;                /* Type GAMBIT de l'élément             */
  int          ndp;                  /* Nombre de noeuds de l'élément        */
  char         str_som_lus[27][9];   /* Sommets lus                          */
  ecs_int_t    som_elt[8];           /* Numeros des noeuds lus de l'élément  */

  ecs_int_t    nbr_som_elt;          /* Nb de noeuds à lire = nbr de sommets */
  ecs_int_t    nbr_elt_ent_max;
  ecs_int_t    nbr_nod_elt_reste;
  ecs_int_t    nbr_nod_ligne;
  ecs_int_t    taille_connect_max;

  ecs_int_t    ind_ent;              /* Indice de l'entité concernée         */

  size_t       ielt;                 /* Indice de boucle sur les éléments    */
  ecs_int_t    inod;                 /* Indice de boucle sur les noeuds      */
  ecs_int_t    inod_ligne;           /* Indice de boucle sur une ligne       */
  ecs_int_t    isom;                 /* Indice de boucle sur les sommets     */
  ecs_int_t    isub;                 /* Indice de boucle sur sous-types      */

  const ecs_loc_gambit_elt_t  *type_elt_gambit;
  const ecs_int_t             *ind_som_elt;     /* Corresp. num. GAMBIT      *
                                                  * / num. locale des         *
                                                  * sommets d'un élément      */

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /*====================================================*/
  /* Initialisations et allocations des tableaux locaux */
  /*====================================================*/

  /* Attention au decalage de `1' !!!         */
  /* On n'alloue pas les tableaux locaux pour */
  /* `ECS_ENTMAIL_DEB = ECS_ENTMAIL_SOM'      */

  for (ind_ent = 0; ind_ent < ECS_N_ENTMAIL; ind_ent++) {

    cpt_elt_ent        [ind_ent] = 0;

    elt_pos_som_ent     [ind_ent] = NULL;
    elt_val_som_ent     [ind_ent] = NULL;
    elt_val_typ_geo_ent [ind_ent] = NULL;
    elt_val_label_ent   [ind_ent] = NULL;

  }

  for (inod = 0; inod < 27; inod++)
    str_som_lus[inod][8] = '\0';

  /*==================================================*/
  /* Boucle de lecture des connectivités des éléments */
  /*==================================================*/

  for (ielt = 0; ielt < nbr_elt; ielt++) {

    ecs_file_gets(chaine, ECS_LOC_LNG_MAX_CHAINE_GAMBIT,
                  fic_maillage, num_ligne);

    /* Lecture des caracteristiques de l'element */
    /*===========================================*/

    /* Format Gambit : I8, 1X, I2, 1X, I2, 1X, 7I8 */

    retour = sscanf(chaine,"%d %d %d",
                    &ne,
                    &ntype,
                    &ndp);

    if (retour != 3)
      ecs_error(__FILE__, __LINE__, 0,
                _("Error reading line %d of file \"%s\"\n"
                  "while decoding element %ld (of %ld)."),
                *num_ligne, ecs_file_get_name(fic_maillage),
                (long)ielt, (long)nbr_elt);

    /* Vérification du type et détermination du sous-type */

    isub = -1;
    type_elt_gambit = NULL;

    if (ntype > 0 && ntype < 8) {
      type_elt_gambit = ecs_loc_gambit_elt_liste_c + ntype - 1;
      for (isub = 0;
           isub < 4 && type_elt_gambit->nbr_som[isub] != ndp;
           isub++);
      if (isub == 4)
        type_elt_gambit = NULL;
    }

    if (type_elt_gambit == NULL)
      ecs_error(__FILE__, __LINE__, 0,
                _("Error reading a GAMBIT mesh file:\n"
                  "at line %d of file \"%s\".\n"
                  "Type identifier <%d> for element <%d> is not recognized."),
                *num_ligne, ecs_file_get_name(fic_maillage),
                (int)ntype, (int)ne);

    /* Lecture des numéros de noeuds constituant l'élément */
    /*=====================================================*/

    /* Format Gambit : 15X, 7I8 (pour chaque ligne) */

    isom     = 0;
    inod     = 0;

    while (inod < ndp) {

      nbr_nod_elt_reste = ndp - inod;
      nbr_nod_ligne     = ECS_MIN(nbr_nod_elt_reste, 7);

      /* Le caractère nul à la fin de chaque chaîne a déjà été positionné */

      for (inod_ligne = 0; inod_ligne < nbr_nod_ligne; inod_ligne++) {
        memcpy(str_som_lus[inod], chaine + 15 + (inod_ligne*8), 8);
        inod += 1;
      }

      /* Passage à une nouvelle ligne si nécessaire */

      if (inod < ndp)
        ecs_file_gets(chaine, ECS_LOC_LNG_MAX_CHAINE_GAMBIT,
                      fic_maillage, num_ligne);

    };

    /* La lecture est terminée, on peut maintenant récuperer les
       numéros des sommets parmi les noeuds */

    nbr_som_elt = ecs_fic_elt_typ_liste_c[type_elt_gambit->ecs_typ].nbr_som;
    ind_som_elt = type_elt_gambit->ind_som[isub];

    for (isom = 0; isom < nbr_som_elt; isom++)
      som_elt[isom] = atoi(str_som_lus[ind_som_elt[isom]]);

    /* Stockage des valeurs lues avant transfert dans la structure `maillage' */
    /*========================================================================*/

    /* Identification de l'entite concernée */
    /*--------------------------------------*/

    ind_ent = ecs_maillage_pre__ret_typ_geo(type_elt_gambit->ecs_typ);

    if (ind_ent != ECS_ENTMAIL_NONE) {

      if (cpt_elt_ent[ind_ent] == 0) {

        /* Si c'est la première fois que l'entité est concernée,
         * on initialise les allocations des tableaux; il ne peut
         * rester plus de (nbr_elt - nombre d'éléments traités)
         * éléments à traiter pour cette entité. */

        nbr_elt_ent_max = nbr_elt - ielt;

        ECS_MALLOC(elt_pos_som_ent[ind_ent], nbr_elt_ent_max + 1, ecs_size_t);

        taille_connect_max = nbr_elt_ent_max*8;
        if (ind_ent == ECS_ENTMAIL_FAC)
          taille_connect_max = nbr_elt_ent_max*4;

        ECS_MALLOC(elt_val_som_ent[ind_ent], taille_connect_max, ecs_int_t);

        ECS_MALLOC(elt_val_typ_geo_ent[ind_ent], nbr_elt_ent_max, ecs_int_t);
        ECS_MALLOC(elt_val_label_ent[ind_ent],   nbr_elt_ent_max, ecs_int_t);

        elt_pos_som_ent[ind_ent][0] = 1;

      }

      /* Connectivité de l'élément par ses numéros de sommets */
      /*------------------------------------------------------*/

      elt_pos_som_ent[ind_ent][cpt_elt_ent[ind_ent] + 1]
        = elt_pos_som_ent[ind_ent][cpt_elt_ent[ind_ent]] + nbr_som_elt;

      for (isom = 0; isom < nbr_som_elt; isom++)
        elt_val_som_ent
          [ind_ent][elt_pos_som_ent[ind_ent][cpt_elt_ent[ind_ent]] - 1 + isom]
          = som_elt[isom];

      /* Type de l'élément lu */
      /*----------------------*/

      elt_val_typ_geo_ent[ind_ent][cpt_elt_ent[ind_ent]]
        = type_elt_gambit->ecs_typ;

      /* Etiquette de l'élément lu */
      /*---------------------------*/

      elt_val_label_ent[ind_ent][cpt_elt_ent[ind_ent]] = ne;

      /* Incrémentation du nombre d'éléments lus */
      /*=========================================*/

      cpt_elt_ent[ind_ent]++;

    }

  } /* Fin de la boucle de lecture des éléments */

  /* Fin de la rubrique */

  ecs_loc_pre_gambit__fin_section(fic_maillage, num_ligne);

  /* Réallocations des tableaux locaux (potentiellement surdimensionnés) */
  /*=====================================================================*/

  for (ind_ent = ECS_ENTMAIL_FAC; ind_ent < ECS_N_ENTMAIL; ind_ent++) {

    if (cpt_elt_ent[ind_ent] != 0) {

      if (cpt_elt_ent[ind_ent] != nbr_elt) {
        ECS_REALLOC(elt_pos_som_ent[ind_ent],
                    cpt_elt_ent[ind_ent] + 1,
                    ecs_size_t);
        ECS_REALLOC(elt_val_typ_geo_ent[ind_ent],
                    cpt_elt_ent[ind_ent],
                    ecs_int_t);
        ECS_REALLOC(elt_val_label_ent[ind_ent],
                    cpt_elt_ent[ind_ent],
                    ecs_int_t);
      }
      ECS_REALLOC(elt_val_som_ent[ind_ent],
                  elt_pos_som_ent[ind_ent][cpt_elt_ent[ind_ent]] - 1,
                  ecs_int_t);

    }
  }
}

/*----------------------------------------------------------------------------
 *  Lecture d'un groupe
 *----------------------------------------------------------------------------*/

static void
ecs_loc_pre_gambit__lit_groupe(ecs_file_t    *fic_maillage,
                               int           *num_ligne,
                               ecs_size_t    *nbr_elt_grp,
                               char         **nom_grp,
                               ecs_int_t    **num_elt_grp)
{
  ecs_int_t  retour;
  ecs_int_t  ind;
  ecs_int_t  ielgp;
  ecs_int_t  ielgp_ligne;
  ecs_int_t  nbr_elgp_reste;
  ecs_int_t  nbr_elgp_ligne;

  /* Variables Gambit lues */

  int  ngp;
  int  nelgp;
  int  mtyp;
  int  nflags;

  char     chaine[ECS_LOC_LNG_MAX_CHAINE_GAMBIT];
  char     ch_num_elt[9];

  const char *mtyp_name[6] = { "Undefined",
                               "Conjugate",
                               "Fluid",
                               "Porous",
                               "Solid",
                               "Deformable" };

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Ligne de description du groupe */

  ecs_file_gets(chaine, ECS_LOC_LNG_MAX_CHAINE_GAMBIT,
                fic_maillage, num_ligne);

  retour = sscanf(chaine, "%*s %d %*s %d %*s %d %*s %d",
                  &ngp, &nelgp, &mtyp, &nflags);

  if (retour != 4)
    ecs_error(__FILE__, __LINE__, 0,
              _("Error reading line %d of file \"%s\"\n"
                "while defining a group."),
              *num_ligne, ecs_file_get_name(fic_maillage));

  if (mtyp < 0 || mtyp > 5)
    mtyp = 0;

  /* Nom du groupe */

  ecs_loc_pre_gambit__raz_chaine(chaine);

  ecs_file_gets(chaine, ECS_LOC_LNG_MAX_CHAINE_GAMBIT,
                fic_maillage, num_ligne);

  chaine[32] = '\0';
  for (ind = 32; ind > 0 && chaine[ind - 1] == ' '; ind--)
    chaine[ind] = '\0';
  for (ind = 0; chaine[ind] != '\0' && chaine[ind] == ' '; ind++);

  /* Stockage du nom de groupe */

  ECS_MALLOC(*nom_grp, strlen(chaine + ind) + 1, char);
  strcpy(*nom_grp, chaine + ind);

  /* Affichage */

  printf(_("\n"
           "  Group %10d: \"%s\"\n"
           "                      type:       %s\n"
           "                      elements:   %d\n"
           "                      indicators: %d\n\n"),
         ngp, chaine + ind, mtyp_name[mtyp], nelgp, nflags);


  /* Saut des indicateurs éventuels (10 / ligne) */

  while (nflags > 0) {
    ecs_file_gets(chaine, ECS_LOC_LNG_MAX_CHAINE_GAMBIT,
                  fic_maillage, num_ligne);
   nflags -= 10;
  }

  /* Initialisations */
  /*=================*/

  ch_num_elt[8] = '\0';

  *nbr_elt_grp = nelgp;

  ECS_MALLOC(*num_elt_grp, nelgp, ecs_int_t);


  /*==========================================*/
  /* Boucle de lecture des numéros de groupes */
  /*==========================================*/

  /* Format Gambit : 10I8 (pour chaque ligne) */

  ielgp = 0;

  while (ielgp < nelgp) {

    /* Lecture d'une nouvelle ligne */

    ecs_file_gets(chaine, ECS_LOC_LNG_MAX_CHAINE_GAMBIT,
                  fic_maillage, num_ligne);

    nbr_elgp_reste = nelgp - ielgp;
    nbr_elgp_ligne= ECS_MIN(nbr_elgp_reste, 10);

    /* Le caractère nul à la fin la chaîne a déjà été positionné */

    for (ielgp_ligne = 0; ielgp_ligne < nbr_elgp_ligne; ielgp_ligne++) {
      memcpy(ch_num_elt, chaine + (ielgp_ligne*8), 8);
      (*num_elt_grp)[ielgp] = atoi(ch_num_elt);
      ielgp += 1;
    }

  };

  /* Fin de la rubrique */

  ecs_loc_pre_gambit__fin_section(fic_maillage, num_ligne);
}

/*----------------------------------------------------------------------------
 *  Lecture d'une condition aux limites
 *----------------------------------------------------------------------------*/

static void
ecs_loc_pre_gambit__lit_cl(ecs_file_t    *fic_maillage,
                           int           *num_ligne,
                           ecs_size_t    *nbr_elt_cl,
                           char         **nom_cl,
                           ecs_int_t    **num_elt_cl,
                           ecs_int_t    **num_fac_cl)
{
  ecs_int_t  retour;
  ecs_int_t  ind;

  int        nbr_codes;
  int        nbr_val_loc;

  /* Variables Gambit lues */

  int  ibcode[5];
  int  nentry;
  int  type;
  int  nvalues;

  char  chaine[ECS_LOC_LNG_MAX_CHAINE_GAMBIT];
  char  ch_num_elt[11];
  char  ch_num_fac[6];
  char  ch_nom_cl[33];
  int   ibcode1_tmp;

  const char *nom_type[3] = { "node",
                              "element/cell",
                              "undefined (incorrect)" };

  const char *nom_ibcode[54] = { "UNSPECIFIED",           /*  0 */
                                 "AXIS",                  /*  1 */
                                 "CONJUGATE",             /*  2 */
                                 "CONVECTION",            /*  3 */
                                 "CYCLIC",                /*  4 */
                                 "DEAD",                  /*  5 */
                                 "ELEMENT_SIDE",          /*  6 */
                                 "ESPECIES",              /*  7 */
                                 "EXHAUST_FAN",           /*  8 */
                                 "FAN",                   /*  9 */
                                 "FREE_SURFACE",          /* 10 */
                                 "GAP",                   /* 11 */
                                 "INFLOW",                /* 12 */
                                 "INLET",                 /* 13 */
                                 "INLET_VENT",            /* 14 */
                                 "INTAKE_FAN",            /* 15 */
                                 "INTERFACE",             /* 16 */
                                 "INTERIOR",              /* 17 */
                                 "INTERNAL",              /* 18 */
                                 "LIVE",                  /* 19 */
                                 "MASS_FLOW_INLET",       /* 20 */
                                 "MELT",                  /* 21 */
                                 "MELT_INTERFACE",        /* 22 */
                                 "MOVING_BOUNDARY",       /* 23 */
                                 "NODE",                  /* 24 */
                                 "OUTFLOW",               /* 25 */
                                 "OUTLET",                /* 26 */
                                 "OUTLET_VENT",           /* 27 */
                                 "PERIODIC",              /* 28 */
                                 "PLOT",                  /* 29 */
                                 "POROUS",                /* 30 */
                                 "POROUS_JUMP",           /* 31 */
                                 "PRESSURE",              /* 32 */
                                 "PRESSURE_FAR_FIELD",    /* 33 */
                                 "PRESSURE_INFLOW",       /* 34 */
                                 "PRESSURE_INLET",        /* 35 */
                                 "PRESSURE_OUTFLOW",      /* 36 */
                                 "PRESSURE_OUTLET",       /* 37 */
                                 "RADIATION",             /* 38 */
                                 "RADIATOR",              /* 39 */
                                 "RECIRCULATION_INLET",   /* 40 */
                                 "RECIRCULATION_OUTLET",  /* 41 */
                                 "SLIP",                  /* 42 */
                                 "SREACTION",             /* 43 */
                                 "SURFACE",               /* 44 */
                                 "SYMMETRY",              /* 45 */
                                 "TRACTION",              /* 46 */
                                 "TRAJECTORY",            /* 47 */
                                 "VELOCITY",              /* 48 */
                                 "VELOCITY_INLET",        /* 49 */
                                 "VENT",                  /* 50 */
                                 "WALL",                  /* 51 */
                                 "SPRING",                /* 52 */
                                 "?" };

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Initialisations */
  /*=================*/

  *nbr_elt_cl = 0;
  *nom_cl     = NULL;
  *num_elt_cl = NULL;

  /* Ligne de description de la condition aux limites */

  ecs_loc_pre_gambit__raz_chaine(chaine);

  ecs_file_gets(chaine, ECS_LOC_LNG_MAX_CHAINE_GAMBIT,
                fic_maillage, num_ligne);

  memcpy(ch_nom_cl, chaine, 32);
  ch_nom_cl[32] = '\0';

  retour = sscanf(chaine + 32, "%d %d %d %d %d %d %d %d",
                  &type, &nentry, &nvalues,
                  ibcode, ibcode+1, ibcode+2, ibcode+3, ibcode+4);

  if (retour < 3)
    ecs_error(__FILE__, __LINE__, 0,
              _("Error reading line %d of file \"%s\"\n"
                "while defining a boundary condition."),
              *num_ligne, ecs_file_get_name(fic_maillage));

  if (type < 0 || type > 1)
    type = 2;

  nbr_codes = retour - 3;

  /* Suppression des blancs en début et fin et stockage du nom */

  for (ind = 32; ind > 0 && ch_nom_cl[ind - 1] == ' '; ind--)
    ch_nom_cl[ind] = '\0';
  for (ind = 0; ch_nom_cl[ind] != '\0' && ch_nom_cl[ind] == ' '; ind++);
  if (ind > 0)
    memmove(ch_nom_cl, ch_nom_cl + ind, 33 - ind);

  /* Affichage */
  /*-----------*/

  printf(_("\n"
           "  Boundary condition: \"%s\"\n"
           "                          type:    %s\n"
           "                          entities: %d\n"
           "                          values:   %d\n"),
         ch_nom_cl, nom_type[type], nentry, nvalues);


  if (nbr_codes > 0) {

    ibcode1_tmp = ibcode[0];
    if (ibcode1_tmp < 0 || ibcode1_tmp > 52)
      ibcode1_tmp = 53;

    printf(_("                          code1:   %d (%s)\n"),
           ibcode[0], nom_ibcode[ibcode[0]]);

  }

  for (ind = 1; ind < nbr_codes; ind++)
    printf(_("                          code%1d:   %d\n"),
           (int)(ind + 1), ibcode[ind]);

  printf("\n");

  /* On ne traite pas les C.L. définies aux sommets */
  /*================================================*/

  if (type != 1) {

    ecs_loc_pre_gambit__saut_section(fic_maillage,
                                     num_ligne);

    return;
  }

  /* Lecture et ajout de la C.L. à la liste */
  /*========================================*/

  *nbr_elt_cl = nentry;

  ECS_MALLOC(*nom_cl, strlen(ch_nom_cl) + 1, char);
  strcpy(*nom_cl, ch_nom_cl);

  ECS_MALLOC(*num_elt_cl, nentry, ecs_int_t);
  ECS_MALLOC(*num_fac_cl, nentry, ecs_int_t);

  /*============================*/
  /* Boucle de lecture des C.L. */
  /*============================*/

  /* Format Gambit : I10, I5, I5 / (4E20.12) */

  ch_num_elt[10] = '\0';
  ch_num_fac[5] = '\0';

  for (ind = 0; ind < nentry; ind++) {

    /* Lecture d'une nouvelle ligne */

    ecs_file_gets(chaine, ECS_LOC_LNG_MAX_CHAINE_GAMBIT,
                  fic_maillage, num_ligne);

    /* Le caractère nul à la fin des chaînes a déjà été positionné */

    memcpy(ch_num_elt, chaine,     10);
    memcpy(ch_num_fac, chaine + 15, 5);

    (*num_elt_cl)[ind] = atoi(ch_num_elt);
    (*num_fac_cl)[ind] = atoi(ch_num_fac);

    /* Saut des valeurs éventuelles (4 / ligne) */

    for (nbr_val_loc = nvalues; nbr_val_loc > 0; nbr_val_loc -= 4)
      ecs_file_gets(chaine, ECS_LOC_LNG_MAX_CHAINE_GAMBIT,
                    fic_maillage, num_ligne);

  };

  /* Fin de la rubrique */

  ecs_loc_pre_gambit__fin_section(fic_maillage, num_ligne);
}

/*----------------------------------------------------------------------------
 *  Transformation de numéros de labels en indice.
 *
 *  On utilise une numérotation de 1 à n, avec un signe positif pour
 *  un indice correspondant à l'entité de plus haut niveau, un signe
 *  négatif pour l'entité inférieure, et 0 pour les labels "non trouvés"
 *  (correspondant probablement à une entité inférieure).
 *----------------------------------------------------------------------------*/

static void
ecs_loc_pre_gambit__label_en_indice(ecs_int_t     nbr_lst,
                                    ecs_size_t    nbr_label_lst[],
                                    ecs_int_t    *val_label_lst[],
                                    size_t        nbr_elt_ent[],
                                    ecs_int_t    *val_label_ent[])
{
  size_t       ind;
  ecs_int_t    ind_lst;
  ecs_int_t    ient;

  ecs_int_t    ient_pcp;

  ecs_int_t    label_max;

  ecs_tab_int_t  tab_label_elt;
  ecs_tab_int_t  tab_label_elt_ord;
  ecs_tab_int_t  tab_ord;

  int sgn;

  const int ient_fac = ECS_ENTMAIL_FAC;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Détermination de l'entité principale et de sa sous-entité */
  /*-----------------------------------------------------------*/

  ient_pcp = ECS_ENTMAIL_FAC;

  if (nbr_elt_ent[ECS_ENTMAIL_CEL] > 0)
    ient_pcp = ECS_ENTMAIL_CEL;

  /* On recherche les labels parmi ceux des éléments */
  /*-------------------------------------------------*/

  for (ient = ient_fac; ient <= ient_pcp; ient++) {

    if (nbr_elt_ent[ient] == 0)
      continue;

    if (ient < ient_pcp)
      sgn = -1;
    else
      sgn = 1;

    tab_label_elt.nbr = nbr_elt_ent[ient];
    tab_label_elt.val = val_label_ent[ient];

    /* On vérifie si la numérotation des labels est dense ou presque */

    label_max = tab_label_elt.val[0];

    for (ind = 1; ind < tab_label_elt.nbr; ind++) {
      if (tab_label_elt.val[ind] > label_max)
        label_max = tab_label_elt.val[ind];
    }

    /* On utilise une recherche par dichotomie */

    tab_ord.nbr = nbr_elt_ent[ient];
    ECS_MALLOC(tab_ord.val, tab_ord.nbr, ecs_int_t);

    tab_label_elt_ord = ecs_tab_int__trie_et_renvoie(tab_label_elt, tab_ord);

    for (ind_lst = 0; ind_lst < nbr_lst; ind_lst++) {

      ecs_tab_int_t  tab_ind;
      ecs_tab_int_t  tab_label;

      tab_label.nbr = nbr_label_lst[ind_lst];
      tab_label.val = val_label_lst[ind_lst];

      if (tab_label.nbr < 1)
        continue;

      tab_ind.nbr = tab_label.nbr;
      ECS_MALLOC(tab_ind.val, tab_ind.nbr, ecs_int_t);

      ecs_tab_int__recherche(tab_label, tab_label_elt_ord, tab_ind);

      /* On transforme les labels en indices */

      for (ind = 0; ind < tab_label.nbr; ind++) {
        if (tab_ind.val[ind] > -1) {
          if (tab_label.val[ind] > 0)
            tab_label.val[ind] = (tab_ord.val[tab_ind.val[ind]] + 1) * sgn;
          else if (ient == ient_pcp)
            tab_label.val[ind] = 0;
        }
      }

      ECS_FREE(tab_ind.val);

    }

    /* Libération mémoire */

    ECS_FREE(tab_label_elt_ord.val);
    ECS_FREE(tab_ord.val);
  }
}

/*----------------------------------------------------------------------------
 *  Construction des éléments surfaciques supplémentaires.
 *  On convertit les références num_elt_cl[] à une entité principale
 *  en des référence à la sous-entité, et on libère les références
 *  num_fac_cl[] à une face de chaque élément de l'entité principale.
 *----------------------------------------------------------------------------*/

static void
ecs_loc_pre_gambit__cree_ent_sub(ecs_int_t    nbr_cl,
                                 ecs_size_t   nbr_elt_cl[],
                                 ecs_int_t   *num_elt_cl[],
                                 ecs_int_t   *num_fac_cl[],
                                 size_t       nbr_elt_ent[],
                                 ecs_size_t  *elt_pos_som_ent[],
                                 ecs_int_t   *elt_val_som_ent[],
                                 ecs_int_t   *elt_val_typ_geo_ent[])
{
  size_t       cpt_elt_sub;
  size_t       taille_connect_sub;

  size_t       ind;
  ecs_int_t    ind_cl;
  ecs_int_t    ind_elt;
  ecs_int_t    ind_fac;
  ecs_int_t    ind_indic;
  ecs_int_t    ind_som;
  ecs_int_t    nbr_som_sub;

  ecs_int_t    type_pcp;
  ecs_int_t    type_sub;

  const ecs_sous_elt_t  *sous_elt;

  ecs_int_t   *val_typ_geo_pcp;
  ecs_int_t   *p_connect_elt;
  ecs_int_t   *p_connect_sub;

  ecs_int_t   *indic_pcp = NULL;

  int  type_gambit[ECS_ELT_TYP_FIN];

  const int ind_ent_cel = ECS_ENTMAIL_CEL;
  const int ind_ent_fac = ECS_ENTMAIL_FAC;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (nbr_cl == 0)
    return;

  if (nbr_elt_ent[ind_ent_cel] == 0)
    return;

  /* Tableau pour correspondance types Enveloppe et GAMBIT */

  for (ind = 0; ind < ECS_ELT_TYP_FIN; ind++)
    type_gambit[ind] = -1;

  for (ind = 0; ind < 7; ind++)
    type_gambit[ecs_loc_gambit_elt_liste_c[ind].ecs_typ] = ind;

  /* Détermination de l'entité principale et de sa sous-entité */
  /*-----------------------------------------------------------*/

  cpt_elt_sub = nbr_elt_ent[ind_ent_fac];

  if (nbr_elt_ent[ind_ent_fac] > 0)
    taille_connect_sub
      = elt_pos_som_ent[ind_ent_fac][nbr_elt_ent[ind_ent_fac]] - 1;
  else
    taille_connect_sub = 0;

  val_typ_geo_pcp = elt_val_typ_geo_ent[ind_ent_cel];

  /* Préparation d'un indicateur sur les cellules (6 faces/cellule max.) */
  /*---------------------------------------------------------------------*/

  ECS_MALLOC(indic_pcp, nbr_elt_ent[ind_ent_cel] * 6, ecs_int_t);

  for (ind = 0; ind < (nbr_elt_ent[ind_ent_cel] * 6); ind++)
    indic_pcp[ind] = 0;

  /* Boucles sur les C.L. pour le comptage */
  /*---------------------------------------*/

  for (ind_cl = 0; ind_cl < nbr_cl; ind_cl++) {

    for (ind = 0; ind < nbr_elt_cl[ind_cl]; ind++) {

      /* On ignore les C.L. sur des sous-sous-entités */

      if (num_elt_cl[ind_cl][ind] > 0) {

        ind_elt = num_elt_cl[ind_cl][ind] - 1;

        /* Ajout "face" si nécessaire */

        ind_indic = ind_elt*6 + num_fac_cl[ind_cl][ind] - 1;

        if (indic_pcp[ind_indic] == 0) {

          ind_fac = num_fac_cl[ind_cl][ind] - 1;
          type_pcp = type_gambit[val_typ_geo_pcp[ind_elt]];
          assert(   ecs_loc_gambit_elt_liste_c[type_pcp].nbr_sselt
                 >= num_fac_cl[ind_cl][ind]);
          sous_elt = &((ecs_loc_gambit_elt_liste_c
                        [type_pcp]).sous_elt[ind_fac]);
          type_sub = sous_elt->elt_typ;

          cpt_elt_sub    += 1;
          taille_connect_sub += (ecs_fic_elt_typ_liste_c[type_sub]).nbr_som;

          indic_pcp[ind_indic] = cpt_elt_sub;

        }
      }
    }
  }

  /* Redimensionnement des tableaux associés à la connectivité */
  /*-----------------------------------------------------------*/

  /* On ne mettra nbr_elt_ent à jour qu'à la fin */

  if (cpt_elt_sub > nbr_elt_ent[ind_ent_fac]) {

    ECS_REALLOC(elt_pos_som_ent[ind_ent_fac], cpt_elt_sub + 1, ecs_size_t);
    ECS_REALLOC(elt_val_som_ent[ind_ent_fac], taille_connect_sub, ecs_int_t);

    if (nbr_elt_ent[ind_ent_fac] == 0)
      elt_pos_som_ent[ind_ent_fac][0] = 1;

  }

  /* Remise des compteurs à leur valeur initiale */

  cpt_elt_sub = nbr_elt_ent[ind_ent_fac];

  if (nbr_elt_ent[ind_ent_fac] > 0)
    taille_connect_sub
      = elt_pos_som_ent[ind_ent_fac][nbr_elt_ent[ind_ent_fac]] - 1;
  else
    taille_connect_sub = 0;

  for (ind = 0; ind < (nbr_elt_ent[ind_ent_cel] * 6); ind++)
    indic_pcp[ind] = 0;

  /* Boucle de construction effective des entités */
  /*-----------------------------------------------*/

  for (ind_cl = 0; ind_cl < nbr_cl; ind_cl++) {

    for (ind = 0; ind < nbr_elt_cl[ind_cl]; ind++) {

      /* On ignore les C.L. sur des sous-sous-entités */

      if (num_elt_cl[ind_cl][ind] > 0) {

        ind_elt = num_elt_cl[ind_cl][ind] - 1;

        /* Ajout "face" si nécessaire */

        ind_fac   = num_fac_cl[ind_cl][ind] - 1;
        ind_indic = ind_elt*6 + ind_fac;

        if (indic_pcp[ind_indic] == 0) {

          type_pcp = type_gambit[val_typ_geo_pcp[ind_elt]];
          assert(  ecs_loc_gambit_elt_liste_c[type_pcp].nbr_sselt > ind_fac);
          sous_elt = &((ecs_loc_gambit_elt_liste_c
                        [type_pcp]).sous_elt[ind_fac]);
          type_sub = sous_elt->elt_typ;

          p_connect_elt =   elt_val_som_ent[ind_ent_cel]
                          + elt_pos_som_ent[ind_ent_cel][ind_elt] - 1;
          p_connect_sub = elt_val_som_ent[ind_ent_fac] + taille_connect_sub;

          nbr_som_sub = (ecs_fic_elt_typ_liste_c[type_sub]).nbr_som;

          for (ind_som = 0; ind_som < nbr_som_sub; ind_som++)
            p_connect_sub[ind_som] = p_connect_elt[sous_elt->som[ind_som]];

          cpt_elt_sub += 1;
          taille_connect_sub += (ecs_fic_elt_typ_liste_c[type_sub]).nbr_som;

          elt_pos_som_ent[ind_ent_fac][cpt_elt_sub] = taille_connect_sub + 1;

          indic_pcp[ind_indic] = cpt_elt_sub;

        }

        /* On convertit la référence à une face de l'entité principale
           en une référence à un élement de la sous-entité */

        num_elt_cl[ind_cl][ind] = - indic_pcp[ind_indic];

      }
      else
        num_elt_cl[ind_cl][ind] = 0;

    }

    /* On n'a plus besoin du tableau num_fac_cl[ind_cl] */

    ECS_FREE(num_fac_cl[ind_cl]);

  }

  nbr_elt_ent[ind_ent_fac] = cpt_elt_sub;

  assert(taille_connect_sub == elt_pos_som_ent[ind_ent_fac][cpt_elt_sub] - 1);

  /* Libération mémoire */
  /*--------------------*/

  for (ind = ECS_ENTMAIL_FAC; ind < ECS_N_ENTMAIL; ind++) {
    if (elt_val_typ_geo_ent[ind] != NULL)
      ECS_FREE(elt_val_typ_geo_ent[ind]);
  }

  ECS_FREE(indic_pcp);
}

/*----------------------------------------------------------------------------
 *  Construction effective des groupes
 *----------------------------------------------------------------------------*/

static void
ecs_loc_pre_gambit__cree_groupes(ecs_maillage_t  *maillage,
                                 ecs_int_t        nbr_grp,
                                 ecs_size_t       nbr_elt_grp[],
                                 char            *nom_grp[],
                                 ecs_int_t       *num_elt_grp[])
{
  ecs_descr_t  *descr_grp; /* Pointeur sur descripteur de table */
  ecs_int_t     num_label; /* Numéro d'indice de l'étiquette */

  size_t       ielgrp;
  size_t       ival;
  int          igrp;
  ecs_int_t    ient;

  bool         bool_aff_grp;
  size_t       nbr_elt_ent[ECS_N_ENTMAIL];
  ecs_int_t    ind_ent;
  ecs_int_t    ind_ent_pcp;

  size_t       ent_cpt_elt[ECS_N_ENTMAIL]; /* Nombre d'elements par entite  */

  /* Stockage avant transfert */
  ecs_int_t  * ent_val_grp[ECS_N_ENTMAIL]; /* Reference /entite et /groupe,
                                              les elts. appart. au groupe  */

  ecs_table_t  *ent_table_grp[ECS_N_ENTMAIL];
  ecs_table_t  *table_grp;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Initialisations */
  /*=================*/

  for (ient = 0; ient < ECS_N_ENTMAIL; ient++) {

    ent_cpt_elt[ient] = 0;
    nbr_elt_ent[ient] = ecs_table__ret_elt_nbr(maillage->table_def[ient]);

    ent_table_grp[ient] = NULL;

  }

  /* Détermination de l'entité principale et de sa sous-entité */
  /*-----------------------------------------------------------*/

  ind_ent_pcp = ECS_ENTMAIL_FAC;

  if (nbr_elt_ent[ECS_ENTMAIL_CEL] > 0)
    ind_ent_pcp = ECS_ENTMAIL_CEL;

  /*===================================*/
  /* Boucle principale sur les groupes */
  /*===================================*/

  for (igrp = 0; igrp < nbr_grp; igrp++) {

    /* On alloue et initialise pour le groupe à traiter */
    /*==================================================*/

    for (ient = ECS_ENTMAIL_FAC; ient < ECS_N_ENTMAIL; ient++) {

      ent_cpt_elt[ient] = 0;

      if (nbr_elt_ent[ient] != 0)
        ECS_MALLOC(ent_val_grp[ient], nbr_elt_ent[ient], ecs_int_t);

      for (ival = 0; ival < nbr_elt_ent[ient]; ival++)
        ent_val_grp[ient][ival] = 0;
    }

    /* traitement d'un groupe */
    /*========================*/

    for (ielgrp = 0; ielgrp < nbr_elt_grp[igrp]; ielgrp++) {

      /* Entité à laquelle correspond un élément du groupe ? */
      /*-----------------------------------------------------*/

      num_label = num_elt_grp[igrp][ielgrp];

      if (num_label > 0)
        ind_ent = ind_ent_pcp;

      else if (num_label < 0) {
        ind_ent = ECS_ENTMAIL_FAC;
        num_label = -num_label;
        assert(ind_ent_pcp > ind_ent);
      }

      else
        ind_ent = ECS_ENTMAIL_NONE;

      if (ind_ent != ECS_ENTMAIL_NONE  && ent_val_grp[ind_ent] != NULL) {

        /* Stockage des valeurs lues avant transfert dans maillage */

        ent_val_grp[ind_ent][num_label - 1] = 1;

        /* Incrémentation du nombre d'objets traités */

        ent_cpt_elt[ind_ent]++;
      }
    }

    /* Boucle de remplissage des entités du maillage */
    /*===============================================*/

    bool_aff_grp = false;

    for (ient = 0; ient < ECS_N_ENTMAIL; ient++) {

      /*--------------------------------------------------------------------*/
      /* S'il y a au moins un élement du groupe de cette entité de maillage */
      /*--------------------------------------------------------------------*/

      if (ent_cpt_elt[ient] != 0) {

        bool_aff_grp = true;

        assert(ent_cpt_elt[ient] <= nbr_elt_ent[ient]);

        /* Création du descripteur de table correspondant au groupe lu */
        /*-------------------------------------------------------------*/

        descr_grp = ecs_descr__cree(ECS_DESCR_IDE_NUL,
                                    nom_grp[igrp]);

        /* Transformation du tableau référencant le groupe en une table */
        /*--------------------------------------------------------------*/

        table_grp = ecs_table__transforme_tableau(nbr_elt_ent[ient],
                                                  ent_val_grp[ient],
                                                  descr_grp);

        if (ent_table_grp[ient] != NULL)
          ecs_table_att__assemble(ent_table_grp[ient],
                                  table_grp);

        else
          ent_table_grp[ient] = table_grp;


      } /* Fin si le nombre d'éléments référencant le groupe n'est pas nul */

    } /* Fin de la boucle sur les entités */


    /* Affichage du bilan des données lues pour les groupes */
    /*======================================================*/

    if (bool_aff_grp == true)
      printf("  %s %d \"%s\"\n",
             _("Group"), igrp + 1, nom_grp[igrp]);

    ecs_maillage_pre__aff_nbr_par_ent(0, ent_cpt_elt, 0);

    /* Incrémentation du compteur sur les groupes */
    /*--------------------------------------------*/

    for (ient = ECS_ENTMAIL_FAC; ient < ECS_N_ENTMAIL; ient++)
      if (nbr_elt_ent[ient] != 0)
        ECS_FREE(ent_val_grp[ient]);

    /* Ré-initialisation des compteurs par entité pour le groupe suivant */
    for (ient = ECS_ENTMAIL_FAC; ient < ECS_N_ENTMAIL; ient++)
      ent_cpt_elt[ient] = 0;

  } /* Fin de la boucle sur les groupes */


  /* Transfert des tables groupe dans les entités de maillage correspondantes */
  /*==========================================================================*/

  for (ient = ECS_ENTMAIL_FAC; ient < ECS_N_ENTMAIL; ient++) {

    assert(maillage->table_att[ient] == NULL);

    if (ent_table_grp[ient] != NULL)
      maillage->table_att[ient] = ent_table_grp[ient];
  }
}

/*============================================================================
 *  Fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Lecture d'un fichier au format GAMBIT neutral
 *   et affectation des donnees dans la structure de maillage
 *
 *  Hypothèses de lecture :
 *   le dataset sur les groupes doit se trouver placer dans le fichier Gambit
 *   après les dataset sur les noeuds et sur les éléments; le dataset
 *   sur noeuds doit lui-même se trouver après celui sur les systèmes de
 *   coordonnées.
 *----------------------------------------------------------------------------*/

ecs_maillage_t *
ecs_pre_gambit__lit_maillage(const char  *nom_fic_maillage)
{
  /* Variables pour le fichier */
  ecs_file_t  *fic_maillage;
  char         chaine[ECS_LOC_LNG_MAX_CHAINE_GAMBIT]; /* Ligne lue */
  int          num_ligne              ; /* Compteur des lignes lues          */

  /* Variables pour la définition des noeuds */
  ecs_coord_t  *som_val_coord;                      /* Coordonnées noeuds    */
  ecs_int_t   *som_val_label;                       /* Labels des noeuds     */

  /* Variables pour la définition des éléments*/
  size_t       cpt_elt_ent        [ECS_N_ENTMAIL]; /* Nombre d'elems/entité */
  ecs_size_t  *elt_pos_som_ent    [ECS_N_ENTMAIL]; /* Positions numéros som */
  ecs_int_t   *elt_val_som_ent    [ECS_N_ENTMAIL]; /* Numéros des sommets   */
  ecs_int_t   *elt_val_typ_geo_ent[ECS_N_ENTMAIL]; /* Types géométriques    */
  ecs_int_t   *elt_val_label_ent  [ECS_N_ENTMAIL]; /* Etiquettes            */

  ecs_int_t    cpt_coul_ent[ECS_N_ENTMAIL];    /* Couleurs (tableaux vides) */
  ecs_int_t   *val_coul_ent[ECS_N_ENTMAIL];
  ecs_size_t  *cpt_elt_coul_ent[ECS_N_ENTMAIL];
  ecs_int_t   *elt_val_color_ent[ECS_N_ENTMAIL];

  /* Autres variables */

  ecs_int_t    ind; /* Indice de boucle */

  bool         bool_elements_lus = false;
  bool         bool_sommets_lus = false;

  int          dim_e   = 3;
  ecs_int_t    nbr_som = 0;
  ecs_int_t    nbr_elt = 0;
  ecs_int_t    nbr_grp = 0;
  ecs_int_t    nbr_cl  = 0;

  ecs_int_t    cpt_grp = 0;
  ecs_int_t    cpt_cl  = 0;

  ecs_size_t  *nbr_elt_grp  = NULL;
  char       **nom_grp      = NULL;
  ecs_int_t  **num_elt_grp  = NULL;

  ecs_size_t  *nbr_elt_cl   = NULL;
  char       **nom_cl       = NULL;
  ecs_int_t  **num_elt_cl   = NULL;
  ecs_int_t  **num_fac_cl   = NULL;

  /* Création d'un maillage initialement vide (valeur de retour) */

  ecs_maillage_t  *maillage = ecs_maillage__cree_nodal();

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Affichage du titre */
  /*====================*/

  printf(_("\n\n"
           "Reading mesh from file in GAMBIT neutral format\n"
           "----------------------\n"));

  printf(_("  Mesh file: %s\n\n\n"), nom_fic_maillage);

  /* Initialisations */
  /*=================*/

  som_val_coord = NULL;
  som_val_label = NULL;

  for (ind = 0; ind < ECS_N_ENTMAIL; ind++) {

    /* Couleurs inutilisées ici, mais les tableaux doivent être définis */
    cpt_coul_ent[ind] = 0;
    val_coul_ent[ind] = NULL;
    cpt_elt_coul_ent[ind] = NULL;
    elt_val_color_ent[ind] = NULL;

  }

  num_ligne = 1;
  dim_e     = 3;

  /* Ouverture du fichier Gambit en lecture */
  /*---------------------------------------*/

  fic_maillage = ecs_file_open(nom_fic_maillage,
                               ECS_FILE_MODE_READ,
                               ECS_FILE_TYPE_TEXT);

  /*================================================*/
  /* Boucle sur les lignes du fichier de maillage : */
  /* tant qu'on n'a pas atteint la fin de fichier   */
  /*================================================*/


  while (ecs_file_gets_try(chaine, ECS_LOC_LNG_MAX_CHAINE_GAMBIT,
                           fic_maillage, &num_ligne) != NULL) {

    /* Décodage selon le type de section de GAMBIT */
    /*=============================================*/

    if (strncmp(chaine, "        CONTROL INFO", 20) == 0) {

      ecs_loc_pre_gambit__lit_entete(fic_maillage,
                                     &num_ligne,
                                     &dim_e,
                                     &nbr_som,
                                     &nbr_elt,
                                     &nbr_grp,
                                     &nbr_cl);

      ECS_MALLOC(nbr_elt_grp, nbr_grp, ecs_size_t);
      ECS_MALLOC(nom_grp, nbr_grp, char *);
      ECS_MALLOC(num_elt_grp, nbr_grp, ecs_int_t *);

      for (ind = 0; ind < nbr_grp; ind++) {
        nom_grp[ind] = NULL;
        num_elt_grp[ind] = NULL;
      }

      ECS_MALLOC(nbr_elt_cl, nbr_cl, ecs_size_t);
      ECS_MALLOC(nom_cl, nbr_cl, char *);
      ECS_MALLOC(num_elt_cl, nbr_cl, ecs_int_t *);
      ECS_MALLOC(num_fac_cl, nbr_cl, ecs_int_t *);

      for (ind = 0; ind < nbr_cl; ind++) {
        nom_cl[ind] = NULL;
        num_elt_cl[ind] = NULL;
        num_fac_cl[ind] = NULL;
      }

    }

    else if (strncmp(chaine, "    APPLICATION DATA", 20) == 0)

      ecs_loc_pre_gambit__saut_section(fic_maillage,
                                       &num_ligne);

    else if (strncmp(chaine, "   NODAL COORDINATES", 20) == 0) {

      /* Lecture des coordonnées des noeuds */

      ECS_MALLOC(som_val_coord, nbr_som * 3, ecs_coord_t);
      ECS_MALLOC(som_val_label, nbr_som, ecs_int_t);


      ecs_loc_pre_gambit__lit_coords(fic_maillage,
                                     dim_e,
                                     nbr_som,
                                     som_val_coord,
                                     som_val_label,
                                     &num_ligne);

      bool_sommets_lus = true;

    }

    else if (strncmp(chaine, "      ELEMENTS/CELLS", 20) == 0) {

      /* Lecture des connectivités des éléments */

      ecs_loc_pre_gambit__lit_elements(fic_maillage,
                                       nbr_elt,
                                       &num_ligne,
                                       cpt_elt_ent,
                                       elt_pos_som_ent,
                                       elt_val_som_ent,
                                       elt_val_typ_geo_ent,
                                       elt_val_label_ent);

      bool_elements_lus = true;

    }

    else if (strncmp(chaine, "       ELEMENT GROUP", 20) == 0) {

      /* Lecture des groupes */

      if (cpt_grp >= nbr_grp)

        ecs_error(__FILE__, __LINE__, 0,
                  _("Error reading a GAMBIT mesh file:\n"
                    "at line %d of file \"%s\".\n"
                    "The number of groups to read is larger\n"
                    "than that defined by the header (%d).n"),
                  num_ligne, ecs_file_get_name(fic_maillage),
                  (int)nbr_grp);


      ecs_loc_pre_gambit__lit_groupe(fic_maillage,
                                     &num_ligne,
                                     &(nbr_elt_grp[cpt_grp]),
                                     &(nom_grp[cpt_grp]),
                                     &(num_elt_grp[cpt_grp]));

      cpt_grp += 1;

    }

    else if (strncmp(chaine, " BOUNDARY CONDITIONS", 20) == 0) {

      /* Lecture des conditions aux limites */

      if (cpt_cl >= nbr_cl)

        ecs_error(__FILE__, __LINE__, 0,
                  _("Error reading a GAMBIT mesh file:\n"
                    "at line %d of file \"%s\".\n"
                    "The number of boundary conditions to read is larger\n"
                    "than that defined by the header (%d).n"),
                  num_ligne, ecs_file_get_name(fic_maillage),
                  (int)nbr_cl);


      ecs_loc_pre_gambit__lit_cl(fic_maillage,
                                 &num_ligne,
                                 &(nbr_elt_cl[cpt_cl]),
                                 &(nom_cl[cpt_cl]),
                                 &(num_elt_cl[cpt_cl]),
                                 &(num_fac_cl[cpt_cl]));

      cpt_cl += 1;

    }

    else if (strncmp(chaine, "   FACE CONNECTIVITY", 20) == 0) {

      ecs_warn();
      printf(_("File \"%s\"\n,"
               "contains a \"FACE CONNECTIVITY\" section\n"
               "indicating non-conforming faces which are not\n"
               "automatically handled by the Preprocessor.\n"
               "-> Use an appropriate joining option.\n"),
             ecs_file_get_name(fic_maillage));

      ecs_loc_pre_gambit__saut_section(fic_maillage, &num_ligne);

    }

    else if (strncmp(chaine, "        TIMESTEPDATA", 20) == 0)

      ecs_loc_pre_gambit__saut_section(fic_maillage, &num_ligne);

    else {

      ecs_warn();
      printf(_("Line %d of file \"%s\"\n,"
               "section: \"%s\" is unknown.\n"),
             num_ligne, ecs_file_get_name(fic_maillage),
             chaine);

      ecs_loc_pre_gambit__saut_section(fic_maillage, &num_ligne);

    }

  }

  /* On vérifie que la fin de fichier a bien été atteinte */

  if (ecs_file_eof(fic_maillage) == 0) {
    ecs_warn();
    printf(_("Line %d of file \"%s\",\n"
             "processing is finished but the end of the file\n"
             "has not been reached.\n\n"),
           num_ligne, ecs_file_get_name(fic_maillage));
  }

  /* On verifie qu'on a bien lu des noeuds et des elements */

  if (bool_sommets_lus == false)
    ecs_error(__FILE__, __LINE__, 0,
              _("Error reading file \"%s\":\n"
                "Section \"NODAL COORDINATES\" was not found."),
              ecs_file_get_name(fic_maillage));

  if (bool_elements_lus == false)
    ecs_error(__FILE__, __LINE__, 0,
              _("Error reading file \"%s\":\n"
                "Section \"ELEMENTS/CELLS\" was not found."),
              ecs_file_get_name(fic_maillage));

  /* Fermeture du fichier de lecture du maillage */
  /*---------------------------------------------*/

  ecs_file_free(fic_maillage);

  /* Construction effective du maillage */
  /*====================================*/

  printf
    (_("  Building mesh structure:\n\n"));

  /* On Convertit les listes à labels en indices */
  /*---------------------------------------------*/

  ecs_loc_pre_gambit__label_en_indice(nbr_grp,
                                      nbr_elt_grp,
                                      num_elt_grp,
                                      cpt_elt_ent,
                                      elt_val_label_ent);

  ecs_loc_pre_gambit__label_en_indice(nbr_cl,
                                      nbr_elt_cl,
                                      num_elt_cl,
                                      cpt_elt_ent,
                                      elt_val_label_ent);

  for (ind = ECS_ENTMAIL_FAC; ind < ECS_N_ENTMAIL; ind++) {
    if (elt_val_label_ent[ind] != NULL)
      ECS_FREE(elt_val_label_ent[ind]);
  }

  /* On peut maintenant construire les éléments de peau associés aux. C.L. */
  /*=======================================================================*/

  ecs_loc_pre_gambit__cree_ent_sub(nbr_cl,
                                   nbr_elt_cl,
                                   num_elt_cl,
                                   num_fac_cl,
                                   cpt_elt_ent,
                                   elt_pos_som_ent,
                                   elt_val_som_ent,
                                   elt_val_typ_geo_ent);


  /* Les conditions aux limites sont maintenant définies comme des
     groupes; on les transforme en groupes définis sur les
     faces au lieu des cellules */

  if (nbr_cl > 0) {

    ECS_REALLOC(nbr_elt_grp, nbr_grp + nbr_cl, ecs_size_t);
    ECS_REALLOC(nom_grp, nbr_grp + nbr_cl, char *);
    ECS_REALLOC(num_elt_grp, nbr_grp + nbr_cl, ecs_int_t *);

    for (ind = 0; ind < nbr_cl; ind++) {
      nbr_elt_grp[nbr_grp + ind] = nbr_elt_cl[ind];
      nom_grp[nbr_grp + ind]     = nom_cl[ind];
      num_elt_grp[nbr_grp + ind] = num_elt_cl[ind];
    }

    nbr_grp += nbr_cl;

  }

  ECS_FREE(nbr_elt_cl);
  ECS_FREE(nom_cl);
  ECS_FREE(num_elt_cl);
  ECS_FREE(num_fac_cl);

  /* Mise à jour des références */
  /*============================*/

  for (ind = ECS_ENTMAIL_FAC; ind < ECS_N_ENTMAIL; ind++) {

    if (cpt_elt_ent[ind] > 0)
      ecs_maillage_pre__label_en_indice
        (nbr_som,
         elt_pos_som_ent[ind][cpt_elt_ent[ind]] - 1,
         som_val_label,
         elt_val_som_ent[ind]);
  }

  ECS_FREE(som_val_label);

  /* On transfère les données dans la structure d'entité de maillage */
  /*=================================================================*/

  ecs_maillage_pre__cree_som(maillage, nbr_som, som_val_coord);

  ecs_maillage_pre__cree_elt(maillage,
                             cpt_elt_ent,
                             elt_pos_som_ent,
                             elt_val_som_ent,
                             NULL,
                             elt_val_color_ent,
                             cpt_coul_ent,
                             val_coul_ent,
                             cpt_elt_coul_ent);

  /* On ajoute les groupes */
  /*-----------------------*/

  ecs_loc_pre_gambit__cree_groupes(maillage,
                                   nbr_grp,
                                   nbr_elt_grp,
                                   nom_grp,
                                   num_elt_grp);

  /* Libérations */
  /*=============*/

  if (nbr_grp > 0) {
    for (ind = 0; ind < nbr_grp; ind++) {
      ECS_FREE(nom_grp[ind]);
      ECS_FREE(num_elt_grp[ind]);
    }
    ECS_FREE(nbr_elt_grp);
    ECS_FREE(nom_grp);
    ECS_FREE(num_elt_grp);
  }

  return maillage;
}

/*----------------------------------------------------------------------------*/

