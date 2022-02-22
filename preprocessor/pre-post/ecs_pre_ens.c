/*============================================================================
 *  Définition de la fonction
 *   de lecture d'un fichier de maillage au format EnSight 6 ou EnSight Gold
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

#include "ecs_pre_ens.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' privés   du  paquetage courant
 *----------------------------------------------------------------------------*/


/*============================================================================
 *                    Déclaration de paramètres et macros
 *============================================================================*/

/* Pour une lecture de 80 caracteres par ligne  */
/* auxquels il faut ajouter le `\n' et le `\0'  */
/* pour l'affectation dans la chaîne réceptrice */
#define ECS_LOC_LNG_MAX_CHAINE_ENS  82


/*============================================================================
 *                  Définition de structures locales
 *============================================================================*/

/* Définition de noeuds */

typedef struct {
  ecs_int_t        nbr_noeuds;  /* Nombre de noeuds */
  ecs_int_t       *id_noeud;    /* Labels des noeuds */
  ecs_coord_t     *coord;       /* Coordonnées (entrelacées) */
  ecs_int_t        num_part;    /* Numéro de part associé */
} ecs_loc_noeuds_ens_t;


/* Définition de connectivités */

typedef struct {
  ecs_int_t        nbr_ele;   /* Nombre d'éléments */
  ecs_elt_typ_t    elt_typ;   /* Type d'élément */
  int32_t         *nbr_n;     /* Nbr. sommets/faces (polygones/polyèdres) */
  int32_t         *nbr_f;     /* Nbr. faces/elem (polyèdres) */
  int32_t         *connect;   /* Connectivité sommets */
  int              num_part;  /* Numéro de part associé */
} ecs_loc_elems_ens_t;


/*============================================================================
 *                              Fonctions privées
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Ouverture et détermination du type du fichier geo
 *----------------------------------------------------------------------------*/

static ecs_file_t *
ecs_loc_pre_ens__ouverture_fic_geo(const char  *nom_fic_geo)
{
  ecs_file_t   *fic;
  ecs_int_t  nbr_char_lus;

  char  chaine[ECS_LOC_LNG_MAX_CHAINE_ENS];

  int32_t  test_endian;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Ouverture du fichier Géométrie et vérification du format */
  /*----------------------------------------------------------*/

  fic = ecs_file_open(nom_fic_geo,
                      ECS_FILE_MODE_READ,
                      ECS_FILE_TYPE_BINARY);

  /* Par défaut, on suppose un fichier binaire C */

  nbr_char_lus = ecs_file_read_try(chaine, sizeof(char), 80, fic);

  /* Si le fichier est trop court pour contenir l'entête binaire, c'est au
     mieux un fichier texte */

  if (nbr_char_lus < 80) {

    ecs_file_set_type(fic, ECS_FILE_TYPE_TEXT);

    ecs_file_rewind(fic);

  }

  /* Pour un fichier binaire Fortran, les 4 premiers octets correspondent à
     la taille de l'enregistrement, et le contenu vient après -> on décale */

  else if (strncmp(chaine + 4,
                   "Fortran Binary",
                   strlen("Fortran Binary")) == 0) {

    ecs_file_set_type(fic, ECS_FILE_TYPE_FORTRAN_BINARY);

    test_endian = *((int32_t *)chaine);

    if (test_endian == 80)
      ecs_file_set_swap_endian(fic, 0);

    else{
      ecs_file_swap_endian(&test_endian, &test_endian, 4, 1);
      if (test_endian == 80)
        ecs_file_set_swap_endian(fic, 1);

      else
        ecs_error(__FILE__, __LINE__, 0,
                  _("EnSight: file \"%s\" seems to be\n"
                    "of Fortran binary type but its header is of the wrong\n"
                    "size or it is of an unknown Fortran binary variant."),
                  ecs_file_get_name(fic));
    }

    /* On se replace après l'entête "Fortran Binary" */

    ecs_file_rewind(fic);
    ecs_file_read(chaine, sizeof(char), 80, fic);

  }

  /* Si le fichier est binaire, on vérifiera plus tard l'aspect
     "big-endian/little-endian" */

  else if (strncmp(chaine, "C Binary", strlen("C Binary")) != 0) {

    ecs_file_set_type(fic, ECS_FILE_TYPE_TEXT);

    ecs_file_rewind(fic);

  }

  return fic;
}

/*----------------------------------------------------------------------------
 *  Lecture d'une chaîne de caractères dans un fichier géométrique EnSight.
 *  Comme les chaînes à lire ne dépassent jamais 80 caractères, on
 *  pourra utiliser un tampon statique.
 *----------------------------------------------------------------------------*/

static void
ecs_loc_pre_ens__lit_chaine(const ecs_file_t  *fic_geo,
                            char   ligne[ECS_LOC_LNG_MAX_CHAINE_ENS],
                            int   *num_ligne)
{
  size_t ind;

  if (ecs_file_get_type(fic_geo) != ECS_FILE_TYPE_TEXT) {

    ecs_file_read(ligne, sizeof(char), 80, fic_geo);
    ligne[80] = '\0';

  }
  else {

    ecs_file_gets(ligne, ECS_LOC_LNG_MAX_CHAINE_ENS, fic_geo, num_ligne);

    for (ind = strlen(ligne) - 1;
         ind > 0 && (ligne[ind] == '\n' || ligne[ind] == '\r');
         ind--);
    ligne[ind + 1] = '\0';

  }

}

/*----------------------------------------------------------------------------
 *  Lecture d'une chaîne de caractères dans un fichier géométrique EnSight.
 *  Comme les chaînes à lire ne dépassent jamais 80 caractères, on
 *  pourra utiliser un tampon statique.
 *  Dans le cas d'un fichier binaire, les chaines à lire sont toujours
 *  exactement de longueur 80 caractères.
 *  Si l'on ne peut lire une chaîne (i.e. fin de fichier), on renvoie NULL;
 *  Sinon, on renvoie un pointeur sur "ligne".
 *----------------------------------------------------------------------------*/

static char *
ecs_loc_pre_ens__lit_chaine_essai(const ecs_file_t  *fic_geo,
                                  char ligne[ECS_LOC_LNG_MAX_CHAINE_ENS] ,
                                  int *num_ligne)
{
  char   *ret;
  size_t  ind;

  if (ecs_file_get_type(fic_geo) != ECS_FILE_TYPE_TEXT) {

    if (ecs_file_read_try(ligne, sizeof(char), 80, fic_geo) == 80) {
      ligne[80] = '\0';
      ret = ligne;
    }
    else
      ret = NULL;

  }
  else {

    ret = ecs_file_gets_try(ligne, ECS_LOC_LNG_MAX_CHAINE_ENS,
                            fic_geo, num_ligne);

    if (ret != NULL) {

      for (ind = strlen(ligne) - 1;
           ind > 0 && (ligne[ind] == '\n' || ligne[ind] == '\r');
           ind--);
      ligne[ind + 1] = '\0';

    }

  }

  return ret;
}

/*----------------------------------------------------------------------------
 *  Lecture d'un tableau d'entiers dans un fichier géométrique EnSight Gold.
 *  Si le tableau val est fourni en argument, on le suppose bien dimensionné
 *  et on le remplit; sinon, si val est à NULL, on alloue un tableau.
 *  On renvoie un pointeur sur le tableau utilisé, qu'il soit fourni en
 *  argument ou alloué ici. Dans ce dernier cas, l'utilisateur aura la
 *  responsabilité de le libérer.
 *----------------------------------------------------------------------------*/

static int32_t *
ecs_loc_pre_ens__lit_int(const ecs_file_t  *fic_geo,
                         ecs_int_t          nbr,
                         int32_t           *val,
                         const int          l_format,
                         int               *num_ligne)
{
  char         ligne[ECS_LOC_LNG_MAX_CHAINE_ENS];
  char         sub[11];
  ecs_int_t    pos_ligne;
  int          ic;
  int          val_lue;

  ecs_int_t    cpt_lus = 0;
  int32_t    * val_loc = val;

  assert(l_format <= 10);

  if (val_loc == NULL)
    ECS_MALLOC(val_loc, nbr, int32_t);


  if (ecs_file_get_type(fic_geo) != ECS_FILE_TYPE_TEXT)
    ecs_file_read(val_loc, sizeof(int32_t), nbr, fic_geo);

  else {

    sub[l_format] = '\0';

    ic = 0;

    while (cpt_lus < nbr) {

      ecs_file_gets(ligne,
                    ECS_LOC_LNG_MAX_CHAINE_ENS,
                    fic_geo,
                    num_ligne);

      /* éventuellement plusieurs valeurs sur la ligne */

      pos_ligne = 0;

      while (   cpt_lus < nbr
             && (   ligne[pos_ligne] != '\0'
                 && ligne[pos_ligne] != '\n'
                 && ligne[pos_ligne] != '\r')) {
        sub[ic] = ligne[pos_ligne];
        ic++;
        pos_ligne++;
        if (ic == l_format) {
          val_lue = atoi(sub);
          val_loc[cpt_lus] = (int32_t)val_lue;
          cpt_lus++;
          ic = 0;
        }
      }

    }

  }

  return val_loc;
}

/*----------------------------------------------------------------------------
 *  Lecture d'un tableau de réels dans un fichier géométrique EnSight.
 *  Si le tableau val est fourni en argument, on le suppose bien dimensionné
 *  et on le remplit ; sinon, si val est à NULL, on alloue un tableau.
 *  On renvoie un pointeur sur le tableau utilisé, qu'il soit fourni en
 *  argument ou alloué ici. Dans ce dernier cas, l'utilisateur aura la
 *  responsabilité de le libérer.
 *----------------------------------------------------------------------------*/

static float *
ecs_loc_pre_ens__lit_float(const ecs_file_t  *fic_geo,
                           ecs_int_t          nbr,
                           float             *val,
                           int               *num_ligne)
{
  char         ligne[ECS_LOC_LNG_MAX_CHAINE_ENS];
  char         sub[13];
  ecs_int_t    pos_ligne;
  int          ic;
  float        val_lue;

  const int        l_format = 12;

  ecs_int_t        cpt_lus = 0;
  float  * val_loc = val;

  if (val_loc == NULL)
    ECS_MALLOC(val_loc, nbr, float);

  if (ecs_file_get_type(fic_geo) != ECS_FILE_TYPE_TEXT)
    ecs_file_read(val_loc, sizeof(float), nbr, fic_geo);

  else {

    sub[l_format] = '\0';

    ic = 0;

    while (cpt_lus < nbr) {

      ecs_file_gets(ligne,
                    ECS_LOC_LNG_MAX_CHAINE_ENS,
                    fic_geo,
                    num_ligne);

      /* éventuellement plusieurs valeurs sur la ligne */

      pos_ligne = 0;

      while (   cpt_lus < nbr
             && (   ligne[pos_ligne] != '\0'
                 && ligne[pos_ligne] != '\n'
                 && ligne[pos_ligne] != '\r')) {
        sub[ic] = ligne[pos_ligne];
        ic++;
        pos_ligne++;
        if (ic == l_format) {
          val_lue = atof(sub);
          val_loc[cpt_lus] = (float)val_lue;
          cpt_lus++;
          ic = 0;
        }
      }

    }

  }

  return val_loc;
}

/*----------------------------------------------------------------------------
 *  Lecture d'un tableau de coordonnées dans un fichier géométrique EnSight 6
 *  (ids et coordonnées mixtes en mode texte : 1 entier et 3 réels par ligne) .
 *  Les tableaux val_id et val_coord sont supposés bien dimensionnés
 *  (val_id étant optionnel).
 *----------------------------------------------------------------------------*/

static void
ecs_loc_pre_ens__lit_tab_coo_6(const ecs_file_t  *fic_geo,
                               ecs_int_t          nbr_lignes,
                               int32_t           *val_id,
                               float             *val_coord,
                               int               *num_ligne)
{
  char         ligne[ECS_LOC_LNG_MAX_CHAINE_ENS];
  char         sub[13];
  ecs_int_t    pos_ligne;
  int          ic;
  int          icoo;

  int          l_format_int = 8;
  int          l_format_real = 12;
  ecs_int_t    cpt_lignes_lues = 0;

  if (ecs_file_get_type(fic_geo) != ECS_FILE_TYPE_TEXT) {

    /* Dans le cas binaire, les tableaux des id (entiers) et coordonnées
       (réels) ne sont pas entrelacés */

    if (val_id != NULL)
      ecs_file_read(val_id, sizeof(int32_t), nbr_lignes, fic_geo);

    ecs_file_read(val_coord, sizeof(float), nbr_lignes * 3, fic_geo);

  }
  else {

    sub[l_format_real] = '\0';

    while (cpt_lignes_lues < nbr_lignes) {

      ecs_file_gets(ligne,
                    ECS_LOC_LNG_MAX_CHAINE_ENS,
                    fic_geo,
                    num_ligne);

      pos_ligne = 0;

      if (val_id != NULL) {
        for (ic = 0; ic < l_format_int; ic++)
          sub[ic] = ligne[pos_ligne++];
        sub[l_format_int] = '\0';
        val_id[cpt_lignes_lues] = (int32_t)(atoi(sub));
      }

      sub[l_format_real] = '\0';
      for (icoo = 0; icoo < 3; icoo++) {
        for (ic = 0; ic < l_format_real; ic++)
          sub[ic] = ligne[pos_ligne++];
        val_coord[cpt_lignes_lues*3 + icoo] = (float)(atof(sub));
      }

      cpt_lignes_lues++;

    }
  }
}

/*----------------------------------------------------------------------------
 *  Lecture d'un tableau de connectivités dans un fichier géométrique
 *  EnSight 6 (ids et connectivités mixtes en mode texte) .
 *  Les tableaux val_id et val_connect sont supposés bien dimensionnés.
 *----------------------------------------------------------------------------*/

static void
ecs_loc_pre_ens__lit_tab_connect(const ecs_file_t  *fic_geo,
                                 ecs_int_t          nbr_lignes,
                                 bool               lit_id,
                                 ecs_int_t          nbr_som_elt,
                                 int32_t           *val_id,
                                 int32_t           *val_connect,
                                 int               *num_ligne)
{
  char         ligne[ECS_LOC_LNG_MAX_CHAINE_ENS];
  char         sub[9];
  ecs_int_t    pos_ligne;
  ecs_int_t    cpt_val_ligne;
  ecs_int_t    cpt_som_ligne;
  int          ic;
  int          val_lue;

  int          l_format = 8;
  ecs_int_t    cpt_lignes_lues = 0;

  if (ecs_file_get_type(fic_geo) != ECS_FILE_TYPE_TEXT) {

    /* Dans le cas binaire, les tableaux des id et les connectivités
       ne sont pas entrelacés */

    if (lit_id == true)
      ecs_file_read(val_id, sizeof(int32_t), nbr_lignes, fic_geo);

    ecs_file_read(val_connect, sizeof(int32_t),
                  nbr_lignes * nbr_som_elt, fic_geo);

  }
  else {

    sub[l_format] = '\0';

    cpt_som_ligne = 0;
    cpt_val_ligne = 0;
    ic = 0;

    while (cpt_lignes_lues < nbr_lignes) { /* Lignes entières, peuvent être
                                              lues via plusieurs appels à
                                              ecs_file_gets si plus
                                              de 80 caractères */

      ecs_file_gets(ligne,
                    ECS_LOC_LNG_MAX_CHAINE_ENS,
                    fic_geo,
                    num_ligne);

      /* éventuellement plusieurs valeurs sur la ligne */

      pos_ligne = 0;

      while (   cpt_lignes_lues < nbr_lignes
             && (   ligne[pos_ligne] != '\0'
                 && ligne[pos_ligne] != '\n'
                 && ligne[pos_ligne] != '\r')) {

        sub[ic] = ligne[pos_ligne];
        ic++;
        pos_ligne++;

        if (ic == l_format) {

          val_lue = atoi(sub);

          if (lit_id == true && cpt_val_ligne == 0)
            val_id[cpt_lignes_lues] = (int32_t)val_lue;
          else {
            val_connect[cpt_lignes_lues*nbr_som_elt + cpt_som_ligne]
              = (int32_t)val_lue;
            cpt_som_ligne++;
          }

          if (cpt_som_ligne == nbr_som_elt) {
            cpt_val_ligne = 0;
            cpt_som_ligne = 0;
            cpt_lignes_lues++;
          }
          else
            cpt_val_ligne++;

          ic = 0;

        }
      }
    }
  }
}

/*----------------------------------------------------------------------------
 *  Lecture des noeuds au format EnSight Gold ;
 *----------------------------------------------------------------------------*/

static void
ecs_loc_pre_ens__lit_noeuds_gold(ecs_file_t             *fic,
                                 ecs_int_t               num_part,
                                 bool                    lire_id_noeud,
                                 bool                    importer_part,
                                 ecs_int_t              *taille_noeuds,
                                 ecs_loc_noeuds_ens_t  **noeuds,
                                 char                    chaine[],
                                 int                    *num_ligne)
{

  ecs_int_t     ind;
  ecs_int_t     ind_coo;
  ecs_int_t     nbr_noeuds;

  int32_t       val_int_lue;

  ecs_loc_noeuds_ens_t  *noeuds_loc;

  int32_t      *id_noeud_loc = NULL;
  ecs_int_t    *id_noeud_tmp = NULL;
  float        *coord_loc    = NULL;
  ecs_coord_t  *coord_tmp    = NULL;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  ecs_loc_pre_ens__lit_chaine(fic, chaine, num_ligne);

  if (strncmp(chaine, "coordinates", strlen("coordinates")) != 0)
    ecs_error(__FILE__, __LINE__, 0,
              _("EnSight: field \"%s\" encountered where\n"
                "\"coordinates\" was expected."), chaine);

  ecs_loc_pre_ens__lit_int(fic, 1, &val_int_lue, 10, num_ligne);

  nbr_noeuds = val_int_lue;

  printf(_("    %10d nodes\n"), (int)nbr_noeuds);

  /* Lecture des labels s'il y a lieu */

  if (lire_id_noeud == true)
    id_noeud_loc = ecs_loc_pre_ens__lit_int(fic, nbr_noeuds,
                                            NULL, 10, num_ligne);

  /* Lecture des coordonnées */

  coord_loc = ecs_loc_pre_ens__lit_float(fic, nbr_noeuds * 3,
                                         NULL, num_ligne);

  if (importer_part == true) {

    ECS_REALLOC(*noeuds, (*taille_noeuds) + 1, ecs_loc_noeuds_ens_t);

    noeuds_loc = (*noeuds) + (*taille_noeuds);

    *taille_noeuds += 1;

    /* Conversion des id noeuds */

    if (id_noeud_loc == NULL || sizeof(int32_t) == sizeof(ecs_int_t))
      id_noeud_tmp = (ecs_int_t *)id_noeud_loc;

    else {
      ECS_MALLOC(id_noeud_tmp, nbr_noeuds, ecs_int_t);
      for (ind = 0; ind < nbr_noeuds; ind++)
        id_noeud_tmp[ind] = id_noeud_loc[ind];
      ECS_FREE(id_noeud_loc);
    }

    /* Entrelacage des coordonnées */

    ECS_MALLOC(coord_tmp, nbr_noeuds*3, ecs_coord_t);

    for (ind = 0; ind < nbr_noeuds; ind++) {
      for (ind_coo = 0; ind_coo < 3; ind_coo++)
        coord_tmp[ind*3 + ind_coo] = coord_loc[nbr_noeuds*ind_coo + ind];
    }

    ECS_FREE(coord_loc);

    noeuds_loc->nbr_noeuds = (ecs_int_t)nbr_noeuds;
    noeuds_loc->id_noeud   = id_noeud_tmp;
    noeuds_loc->coord      = coord_tmp;
    noeuds_loc->num_part   = num_part;

  }
  else {

    ECS_FREE(coord_loc);
    ECS_FREE(id_noeud_loc);

  }
}

/*----------------------------------------------------------------------------
 *  Lecture des noeuds au format EnSight 6 ;
 *----------------------------------------------------------------------------*/

static void
ecs_loc_pre_ens__lit_noeuds_6(ecs_file_t             *fic,
                              bool                    lire_id_noeud,
                              ecs_int_t              *taille_noeuds,
                              ecs_loc_noeuds_ens_t  **noeuds,
                              char                    chaine[],
                              int                    *num_ligne)
{
  ecs_int_t    ind;
  ecs_int_t    nbr_noeuds;

  int32_t      val_int_lue;

  int32_t       *id_noeud_loc = NULL;
  ecs_int_t     *id_noeud_tmp = NULL;
  float         *coord_loc    = NULL;
  ecs_coord_t   *coord_tmp    = NULL;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  ecs_loc_pre_ens__lit_chaine(fic, chaine, num_ligne);

  if (strncmp(chaine, "coordinates", strlen("coordinates")) != 0)
    ecs_error(__FILE__, __LINE__, 0,
              _("EnSight: field \"%s\" encountered where\n"
                "\"coordinates\" was expected."), chaine);

  ecs_loc_pre_ens__lit_int(fic, 1, &val_int_lue, 8, num_ligne);

  nbr_noeuds = val_int_lue;

  if (lire_id_noeud == true)
    ECS_MALLOC(id_noeud_loc, nbr_noeuds,  int32_t);
  else
    id_noeud_loc = NULL;

  ECS_MALLOC(coord_loc, nbr_noeuds*3,  float);

  ecs_loc_pre_ens__lit_tab_coo_6(fic,
                                 nbr_noeuds ,
                                 id_noeud_loc ,
                                 coord_loc ,
                                 num_ligne);

  if (coord_loc != NULL) {
    *taille_noeuds = 1;
    ECS_MALLOC(*noeuds, (*taille_noeuds),  ecs_loc_noeuds_ens_t);
  }
  else
    ecs_error(__FILE__, __LINE__, 0,
              _("EnSight: no node defined in file\n\"%s\"."),
              ecs_file_get_name(fic));

  /* Conversion des id noeuds */

  if (id_noeud_loc == NULL || sizeof(int32_t) == sizeof(ecs_int_t))
    id_noeud_tmp = (ecs_int_t *)id_noeud_loc;

  else {
    ECS_MALLOC(id_noeud_tmp, nbr_noeuds, ecs_int_t);
    for (ind = 0; ind < nbr_noeuds; ind++)
      id_noeud_tmp[ind] = id_noeud_loc[ind];
    ECS_FREE(id_noeud_loc);
  }

  /* Conversion des coordonnées */

  ECS_MALLOC(coord_tmp, nbr_noeuds*3, ecs_coord_t);

  for (ind = 0; ind < 3*nbr_noeuds; ind++) {
    coord_tmp[ind] = coord_loc[ind];
  }

  ECS_FREE(coord_loc);

  (*noeuds)->nbr_noeuds = (ecs_int_t)nbr_noeuds;
  (*noeuds)->id_noeud   = id_noeud_tmp;
  (*noeuds)->coord      = coord_tmp;
  (*noeuds)->num_part   =  -1;
}

/*----------------------------------------------------------------------------
 * Retourne le type de l'élément ainsi que son nombre de sommet
 *----------------------------------------------------------------------------*/

static void
ecs_loc_pre_ens__type_elem(ecs_file_t     *fic,
                           const char      chaine[],
                           int            *nbr_som_elt,
                           ecs_elt_typ_t  *elt_typ)
{
  /*----------------------------*/
  /* Décodage du type d'élément */
  /*----------------------------*/

  if (strncmp(chaine, "g_", strlen("g_")) == 0)
    ecs_error(__FILE__, __LINE__, 0,
              _("EnSight: error reading file \"%s\".\n"
                "The current \"part\" contains ghost cells,\n"
                "not currently handled."),
              ecs_file_get_name(fic));

  else if (strncmp(chaine, "bar", strlen("bar")) == 0) {
    *elt_typ         = ECS_ELT_TYP_NUL;
    *nbr_som_elt     = atoi(chaine+strlen("bar"));
  }
  else if (strncmp(chaine, "tria", strlen("tria")) == 0) {
    *elt_typ         = ECS_ELT_TYP_FAC_TRIA;
    *nbr_som_elt     = atoi(chaine+strlen("tria"));
  }
  else if (strncmp(chaine, "quad", strlen("quad")) == 0) {
    *elt_typ         = ECS_ELT_TYP_FAC_QUAD;
    *nbr_som_elt     = atoi(chaine+strlen("quad"));
  }
  else if (strncmp(chaine, "tetra", strlen("tetra")) == 0) {
    *elt_typ         = ECS_ELT_TYP_CEL_TETRA;
    *nbr_som_elt     = atoi(chaine+strlen("tetra"));
  }
  else if (strncmp(chaine, "pyramid", strlen("pyramid")) == 0) {
    *elt_typ         = ECS_ELT_TYP_CEL_PYRAM;
    *nbr_som_elt     = atoi(chaine+strlen("pyramid"));
  }
  else if (strncmp(chaine, "penta", strlen("penta")) == 0) {
    *elt_typ         = ECS_ELT_TYP_CEL_PRISM;
    *nbr_som_elt     = atoi(chaine+strlen("penta"));
  }
  else if (strncmp(chaine, "hexa", strlen("hexa")) == 0) {
    *elt_typ         = ECS_ELT_TYP_CEL_HEXA;
    *nbr_som_elt     = atoi(chaine+strlen("hexa"));
  }
  else if (strncmp(chaine, "nsided", strlen("nsided")) == 0) {
    *elt_typ         = ECS_ELT_TYP_FAC_POLY;
    *nbr_som_elt     = 0;
  }
  else if (strncmp(chaine, "nfaced", strlen("nfaced")) == 0) {
    *elt_typ         = ECS_ELT_TYP_CEL_POLY;
    *nbr_som_elt     = 0;
  }
  else {
    ecs_error(__FILE__, __LINE__, 0,
              _("EnSight: error reading file \"%s\".\n"
                "The current \"part\" contains a section of type:\n"
                "\"%s\"\n"
                "not currently handled."),
              ecs_file_get_name(fic), chaine);
  }
}

/*----------------------------------------------------------------------------
 * Transformation des éléments en éléments linéaires.
 * Cette fonction renvoie un pointeur sur le tableau de connectivité réalloué.
 *----------------------------------------------------------------------------*/

static int32_t *
ecs_loc_pre_ens__ele_lin(int32_t     *connect_loc,
                         ecs_int_t    nbr_ele,
                         ecs_int_t    taille_ele,
                         ecs_int_t    taille_ele_lin)
{
  ecs_int_t    ielt;
  ecs_int_t    iloc;
  ecs_int_t    ipos_lin;
  ecs_int_t    ipos_ens;

  ipos_lin = 0;

  for (ielt = 0; ielt < nbr_ele; ielt++) {

    ipos_ens = ielt * taille_ele;

    for (iloc = 0; iloc < taille_ele_lin; iloc++)
      connect_loc[ipos_lin++] = connect_loc[ipos_ens++];

  }

  ECS_REALLOC(connect_loc, ipos_lin, int32_t);

  return connect_loc;
}

/*----------------------------------------------------------------------------
 * Lecture des éléments format EnSight Gold; renvoie le nombre d'éléments lus,
 * ou -1 si l'on a affaire à une "part"
 *----------------------------------------------------------------------------*/

static ecs_int_t
ecs_loc_pre_ens__lit_elem_gold(ecs_file_t                   *fic,
                               ecs_int_t                     num_part,
                               bool                          lire_id_elem,
                               bool                          importer_part,
                               const ecs_loc_noeuds_ens_t   *noeuds,
                               ecs_int_t                    *taille_elems,
                               ecs_loc_elems_ens_t         **elems,
                               char                          chaine[],
                               int                          *num_ligne)
{
  int            nbr_som_elt;

  ecs_int_t      ind;
  ecs_int_t      taille_connect;
  ecs_int_t      taille_lect;

  ecs_elt_typ_t   elt_typ;

  int32_t     * id_elem;

  ecs_loc_elems_ens_t  *elems_loc;

  ecs_int_t  nbr_som_elt_lin = 0;
  int32_t    nbr_elt_loc = 0;
  int32_t  * nbr_n_loc = NULL;
  int32_t  * nbr_f_loc = NULL;
  int32_t  * connect_loc = NULL;

  const ecs_int_t  l_fmt_int = 10;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* On s'assure qu'on n'a pas affaire à une "part" */

  if (strncmp(chaine, "part", strlen("part")) == 0)
    return -1;

  /*----------------------------*/
  /* Décodage du type d'élément */
  /*----------------------------*/

  ecs_loc_pre_ens__type_elem(fic, chaine, &nbr_som_elt, &elt_typ);

  /* Lecture du nombre d'éléments */
  /*------------------------------*/

  ecs_loc_pre_ens__lit_int(fic, 1, &nbr_elt_loc, l_fmt_int, num_ligne);

  printf(_("    %10d %-10s elements\n"), (int)nbr_elt_loc, chaine);

  /*------------------------------------------------*/
  /* Lecture éventuelle des labels des éléments     */
  /*------------------------------------------------*/

  if (lire_id_elem == true) {
    id_elem = ecs_loc_pre_ens__lit_int(fic,
                                       (ecs_int_t)nbr_elt_loc,
                                       NULL,
                                       l_fmt_int,
                                       num_ligne);
    ECS_FREE(id_elem);
  }

  /*------------------------------------------------*/
  /* Lecture de la connectivité nodale des elements */
  /*------------------------------------------------*/

  if (elt_typ == ECS_ELT_TYP_FAC_POLY) {

    nbr_n_loc = ecs_loc_pre_ens__lit_int(fic,
                                         (ecs_int_t)nbr_elt_loc,
                                         NULL,
                                         l_fmt_int,
                                         num_ligne);
    taille_lect = 0;
    for (ind = 0; ind < (ecs_int_t)nbr_elt_loc; ind++)
      taille_lect += nbr_n_loc[ind];

    taille_connect = taille_lect;

  }
  else if (elt_typ == ECS_ELT_TYP_CEL_POLY) {

    ecs_int_t taille_int = 0;

    nbr_f_loc = ecs_loc_pre_ens__lit_int(fic,
                                         (ecs_int_t)nbr_elt_loc,
                                         NULL,
                                         l_fmt_int,
                                         num_ligne);
    taille_lect = 0;

    for (ind = 0; ind < (ecs_int_t)nbr_elt_loc; ind++)
      taille_int += nbr_f_loc[ind];

    nbr_n_loc = ecs_loc_pre_ens__lit_int(fic,
                                         taille_int,
                                         NULL,
                                         l_fmt_int,
                                         num_ligne);

    for (ind = 0; ind < taille_int; ind++)
      taille_lect += nbr_n_loc[ind];

    taille_connect = taille_lect;

  }
  else  {
    nbr_som_elt_lin = ecs_fic_elt_typ_liste_c[elt_typ].nbr_som;
    taille_lect    = nbr_elt_loc * nbr_som_elt;
    taille_connect = nbr_elt_loc * nbr_som_elt_lin;
  }
  connect_loc = ecs_loc_pre_ens__lit_int(fic,
                                         taille_lect,
                                         NULL,
                                         l_fmt_int,
                                         num_ligne);

  if (importer_part == true && elt_typ != ECS_ELT_TYP_NUL) {

    ECS_REALLOC(*elems, (*taille_elems) + 1, ecs_loc_elems_ens_t);

    elems_loc = (*elems) + (*taille_elems);

    *taille_elems += 1;

    /* Suppression références noeuds non sommets éventuels */

    if (taille_lect > taille_connect) {
      connect_loc = ecs_loc_pre_ens__ele_lin(connect_loc,
                                             (ecs_int_t)nbr_elt_loc,
                                             nbr_som_elt,
                                             nbr_som_elt_lin);

    }

    /* Application des labels des sommets */

    if (noeuds->id_noeud != NULL) {

      for (ind = 0; ind < taille_connect; ind++){
        assert(   (connect_loc[ind]-1) >= 0
               && (connect_loc[ind]-1) < noeuds->nbr_noeuds);
        connect_loc[ind] = noeuds->id_noeud[connect_loc[ind]-1];
      }
    }

    elems_loc->nbr_ele  = nbr_elt_loc;
    elems_loc->elt_typ  = elt_typ;
    elems_loc->nbr_n    = nbr_n_loc;
    elems_loc->nbr_f    = nbr_f_loc;
    elems_loc->connect  = connect_loc;
    elems_loc->num_part = num_part;

  }
  else {

    if (nbr_n_loc != NULL)
      ECS_FREE(nbr_n_loc);

    if (nbr_f_loc != NULL)
      ECS_FREE(nbr_f_loc);

    ECS_FREE(connect_loc);

  }

  return nbr_elt_loc;

}

/*----------------------------------------------------------------------------
 * Lecture des éléments format EnSight 6; renvoie le nombre d'éléments lus,
 * ou -1 si l'on a affaire à une "part"
 *----------------------------------------------------------------------------*/

static ecs_int_t
ecs_loc_pre_ens__lit_elem_6(ecs_file_t            *fic,
                            int                    num_part,
                            bool                   lire_id_elem,
                            bool                   importer_part,
                            ecs_int_t             *taille_elems,
                            ecs_loc_elems_ens_t  **elems,
                            char                   chaine[],
                            int                   *num_ligne)
{
  int            nbr_som_elt;

  ecs_int_t      nbr_som_elt_lin;
  ecs_int_t      taille_connect;
  ecs_int_t      taille_lect;

  ecs_elt_typ_t   elt_typ;

  int32_t  * id_elem;

  ecs_loc_elems_ens_t  *elems_loc;

  int32_t    nbr_elt_loc = 0;
  int32_t  * nbr_n_loc = NULL;
  int32_t  * nbr_f_loc = NULL;
  int32_t  * connect_loc = NULL;

  const ecs_int_t  l_fmt_int = 8;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* On s'assure qu'on n'a pas affaire à une "part" */

  if (strncmp(chaine, "part", strlen("part")) == 0)
    return -1;

  /*----------------------------*/
  /* Décodage du type d'élément */
  /*----------------------------*/

  ecs_loc_pre_ens__type_elem(fic, chaine, &nbr_som_elt, &elt_typ);

  /* Lecture du nombre d'éléments */
  /*------------------------------*/

  ecs_loc_pre_ens__lit_int(fic, 1, &nbr_elt_loc, l_fmt_int, num_ligne);

  printf(_("    %10d %-10s elements\n"), (int)nbr_elt_loc, chaine);

  /*------------------------------------------------*/
  /* Lecture éventuelle des labels des éléments     */
  /* et Lecture de la connectivité                  */
  /*------------------------------------------------*/

  ECS_MALLOC(id_elem, nbr_elt_loc, int32_t);
  ECS_MALLOC(connect_loc, nbr_elt_loc * nbr_som_elt, int32_t);

  ecs_loc_pre_ens__lit_tab_connect(fic ,
                                   nbr_elt_loc  ,
                                   lire_id_elem ,
                                   nbr_som_elt ,
                                   id_elem ,
                                   connect_loc ,
                                   num_ligne);
  ECS_FREE(id_elem);

  nbr_som_elt_lin = ecs_fic_elt_typ_liste_c[elt_typ].nbr_som;
  taille_lect    = nbr_elt_loc * nbr_som_elt;
  taille_connect = nbr_elt_loc * nbr_som_elt_lin;

  if (importer_part == true && elt_typ != ECS_ELT_TYP_NUL) {

    ECS_REALLOC(*elems, (*taille_elems) + 1, ecs_loc_elems_ens_t);

    elems_loc = (*elems) + (*taille_elems);

    *taille_elems += 1;

    /* Suppression références noeuds non sommets éventuels */

    if (taille_lect > taille_connect) {
      connect_loc = ecs_loc_pre_ens__ele_lin(connect_loc,
                                             (ecs_int_t)nbr_elt_loc,
                                             nbr_som_elt,
                                             nbr_som_elt_lin);

    }

    elems_loc->nbr_ele  = nbr_elt_loc;
    elems_loc->elt_typ  = elt_typ;
    elems_loc->nbr_n    = nbr_n_loc;
    elems_loc->nbr_f    = nbr_f_loc;
    elems_loc->connect  = connect_loc;
    elems_loc->num_part = num_part;

  }
  else {

    ECS_FREE(connect_loc);

  }
  return nbr_elt_loc;
}

/*----------------------------------------------------------------------------
 *  Vérification qu'un groupe de noeuds est bien défini avec des
 *  labels croissants, et tri le cas échéant.
 *----------------------------------------------------------------------------*/

static void
ecs_loc_pre_ens__trie_noeuds(ecs_loc_noeuds_ens_t  *noeuds)
{
  ecs_int_t         ind;
  ecs_int_t         ind_coo;

  ecs_tab_int_t     id_noeud_loc;
  ecs_tab_int_t     id_trie;
  ecs_tab_int_t     renum;

  bool              a_trier = false;

  ecs_coord_t      *coord_tmp    = NULL;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (noeuds->id_noeud == NULL)
    return;

  a_trier = false;

  for (ind = 0; ind < noeuds->nbr_noeuds - 1; ind++) {

    if (noeuds->id_noeud[ind + 1] <= noeuds->id_noeud[ind]) {

      a_trier = true;
      break;

    }

  }

  /* Si l'on a des sommets à trier */

  if (a_trier == true) {

    id_noeud_loc.nbr = noeuds->nbr_noeuds;
    id_noeud_loc.val = noeuds->id_noeud;

    /* Tri effectif */

    ECS_MALLOC(renum.val, noeuds->nbr_noeuds, ecs_int_t);
    renum.nbr = noeuds->nbr_noeuds;

    id_trie = ecs_tab_int__trie_et_renvoie(id_noeud_loc,
                                           renum);

    ECS_FREE(noeuds->id_noeud);

    noeuds->id_noeud = id_trie.val;

    ECS_MALLOC(coord_tmp, noeuds->nbr_noeuds*3, ecs_coord_t);

    for (ind = 0; ind < noeuds->nbr_noeuds; ind++) {
      for (ind_coo = 0; ind_coo < 3; ind_coo++)
        coord_tmp[ind*3 + ind_coo]
          = noeuds->coord[(renum.val[ind])*3 + ind_coo];

    }

    ECS_FREE(renum.val);

    ECS_FREE(noeuds->coord);

    noeuds->coord = coord_tmp;
 }
}

/*----------------------------------------------------------------------------
 *  Fusion de deux blocs de définitions de noeuds ; le premier ensemble
 *  (noeuds_ref) est étendu au besoin, le deuxième est supprimé (i.e.
 *  ses tableaux sont libérés, et il n'est plus valide à l'issue de
 *  cette fonction).
 *----------------------------------------------------------------------------*/

static void
ecs_loc_pre_ens__fusion_noeuds(ecs_loc_noeuds_ens_t  *noeuds_ref,
                               ecs_loc_noeuds_ens_t  *noeuds_add)
{
  ecs_int_t          cpt_add;
  ecs_int_t          ind_add;
  ecs_int_t          ind_ref;
  ecs_int_t          ind_tot;
  ecs_int_t          ind_coo;
  ecs_int_t          nbr_noeuds_tot;

  ecs_int_t        * id_new    = NULL;
  ecs_coord_t      * coord_new = NULL;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (noeuds_ref->id_noeud == NULL || noeuds_ref->id_noeud == NULL)
    ecs_error(__FILE__, __LINE__, 0,
              _("EnSight: error merging node groups.\n"
                "At least one of the groups has no \"id\", and they\n"
                "should not be merged (internal logic problem)."));


  ind_ref = 0;
  cpt_add = 0;

  for (ind_add = 0; ind_add < noeuds_add->nbr_noeuds; ind_add++) {

    while (   ind_ref < noeuds_ref->nbr_noeuds
           && noeuds_ref->id_noeud[ind_ref] < noeuds_add->id_noeud[ind_add])
      ind_ref++;

    if (   ind_ref == noeuds_ref->nbr_noeuds
        || (ind_ref < noeuds_ref->nbr_noeuds
            && (   noeuds_ref->id_noeud[ind_ref]
                != noeuds_add->id_noeud[ind_add])))
      cpt_add++;

  }

  /* Si l'on a des sommets à fusionner */

  if (cpt_add > 0) {

    nbr_noeuds_tot = noeuds_ref->nbr_noeuds + cpt_add;

    ECS_MALLOC(id_new,    nbr_noeuds_tot,   ecs_int_t);
    ECS_MALLOC(coord_new, nbr_noeuds_tot*3, ecs_coord_t);

    ind_ref = 0;
    ind_tot = 0;
    cpt_add = 0;

    for (ind_add = 0; ind_add < noeuds_add->nbr_noeuds; ind_add++) {

      while (   ind_ref < noeuds_ref->nbr_noeuds
             && noeuds_ref->id_noeud[ind_ref] < noeuds_add->id_noeud[ind_add]) {
        id_new[ind_tot] = noeuds_ref->id_noeud[ind_ref];
        for (ind_coo = 0; ind_coo < 3; ind_coo++)
          coord_new[ind_tot*3 + ind_coo]
            = noeuds_ref->coord[ind_ref*3 + ind_coo];
        ind_tot++;
        ind_ref++;
      }

      if (   ind_ref == noeuds_ref->nbr_noeuds
          || (ind_ref < noeuds_ref->nbr_noeuds
              && (   noeuds_ref->id_noeud[ind_ref]
                  != noeuds_add->id_noeud[ind_add]))) {
        id_new[ind_tot] = noeuds_add->id_noeud[ind_add];
        for (ind_coo = 0; ind_coo < 3; ind_coo++)
          coord_new[ind_tot*3 + ind_coo]
            = noeuds_add->coord[ind_add*3 + ind_coo];
        ind_tot++;
      }

    }

    /* Il peut rester des éléments dans la liste de référence */

    while (   ind_ref < noeuds_ref->nbr_noeuds ) {
      id_new[ind_tot] = noeuds_ref->id_noeud[ind_ref];
      for (ind_coo = 0; ind_coo < 3; ind_coo++)
        coord_new[ind_tot*3 + ind_coo]
          = noeuds_ref->coord[ind_ref*3 + ind_coo];
      ind_tot++;
      ind_ref++;
    }

    ECS_FREE(noeuds_ref->id_noeud);
    ECS_FREE(noeuds_ref->coord);

    /* On remplace les anciennes valeurs de référence par les
       valeurs fusionnées */

    noeuds_ref->nbr_noeuds = ind_tot;
    noeuds_ref->id_noeud   = id_new;
    noeuds_ref->coord      = coord_new;

  }

  /* On libère les valeurs à fusionner */

  ECS_FREE(noeuds_add->id_noeud);
  ECS_FREE(noeuds_add->coord);
}

/*----------------------------------------------------------------------------
 *  Concaténation des éléments pour préparer le transfert dans la structure
 *   de maillage ; la liste des élémens liste_elems est libérée, et remplacée
 *   par la structure renvoyée.
 *
 *  Les références aux ids des noeuds sont transformés en indices, et
 *   les ids des noeuds ensuite supprimés.
 *----------------------------------------------------------------------------*/

static void
ecs_loc_pre_ens__concat_elems(ecs_maillage_t         *maillage,
                              ecs_int_t               taille_liste_elems,
                              ecs_loc_noeuds_ens_t   *liste_noeuds,
                              ecs_loc_elems_ens_t   **liste_elems)
{
  ecs_entmail_t    entmail_e;

  ecs_int_t      ient;
  ecs_int_t      ielt;
  ecs_int_t      ielts_loc;
  ecs_int_t      ifac;
  ecs_int_t      isom;
  ecs_int_t      nbr_som_elt;
  ecs_int_t      nbr_som_fac;
  ecs_int_t      pos_elt;

  ecs_loc_elems_ens_t  * elems_loc;
  ecs_elt_typ_t          typ_geo;

  /* Déclarations des variables de stockage        */
  /* avant transfert dans la structure du maillage */
  /*-----------------------------------------------*/

  size_t       cpt_elt_ent        [ECS_N_ENTMAIL]; /* Nbr. elts par entité */
  ecs_int_t    cpt_coul_ent       [ECS_N_ENTMAIL]; /* Compteur de couleurs */
  ecs_int_t    cpt_val_som_ent    [ECS_N_ENTMAIL]; /* Taille connect.      */
  ecs_int_t    icoul              [ECS_N_ENTMAIL]; /* Couleur en cours     */
  ecs_int_t  * val_coul_ent       [ECS_N_ENTMAIL]; /* Tableau des couleurs */
  ecs_size_t * cpt_elt_coul_ent   [ECS_N_ENTMAIL];
  ecs_size_t * elt_pos_som_ent    [ECS_N_ENTMAIL]; /* Positions sommets    */
  ecs_int_t  * elt_val_som_ent    [ECS_N_ENTMAIL]; /* Numéros des sommets  */
  ecs_int_t  * elt_val_color_ent  [ECS_N_ENTMAIL]; /* Couleurs éléments    */

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /*====================================================*/
  /* Initialisations et allocations des tableaux locaux */
  /*====================================================*/

  /* Attention au decalage de `1' !!!         */
  /* On n'alloue pas les tableaux locaux pour */
  /* `ECS_ENTMAIL_DEB = ECS_ENTMAIL_SOM'      */

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

  /* Vérification que la liste des éléments n'est pas vide */

  assert(taille_liste_elems != 0);
  assert(liste_elems != NULL);

  /* Boucle sur la liste des éléments pour le dimensionnement */
  /*----------------------------------------------------------*/

  for (ielts_loc = 0; ielts_loc < taille_liste_elems; ielts_loc++) {

    elems_loc = (*liste_elems) + ielts_loc;

    /* Identification de l'entité concernée */

    typ_geo = elems_loc->elt_typ;
    entmail_e = ecs_maillage_pre__ret_typ_geo(elems_loc->elt_typ);

    /* Nombre d'éléments à ajouter */

    cpt_elt_ent[entmail_e] += elems_loc->nbr_ele;

    /* Traitement des éléments "classiques" */

    if (typ_geo != ECS_ELT_TYP_CEL_POLY && typ_geo != ECS_ELT_TYP_FAC_POLY) {

      nbr_som_elt = ecs_fic_elt_typ_liste_c[typ_geo].nbr_som;

      cpt_val_som_ent[entmail_e] += (elems_loc->nbr_ele * nbr_som_elt);

    }

    /* Traitement des Polygones */

    else if (typ_geo == ECS_ELT_TYP_FAC_POLY) {

      for (ielt = 0; ielt < elems_loc->nbr_ele; ielt++)

        cpt_val_som_ent[entmail_e] += (elems_loc->nbr_n[ielt]);

    }

    /* Traitement des Polyèdres */

    else if (typ_geo == ECS_ELT_TYP_CEL_POLY) {

      ecs_int_t  cpt_fac_loc = 0;

      for (ielt = 0; ielt < elems_loc->nbr_ele; ielt++) {

        for (ifac = 0; ifac < elems_loc->nbr_f[ielt]; ifac++) {

          /* Le premier noeud de chaque face est répété en queue pour la
             détection de fin de face -> + 1 pour la taille */

          cpt_val_som_ent[entmail_e] += elems_loc->nbr_n[cpt_fac_loc] + 1;

          cpt_fac_loc++;

        }

      }

    } /* Fin du traitement selon le type d'entité */

  }

  /* Allocations et initialisation */
  /*-------------------------------*/

  for (ient = ECS_ENTMAIL_FAC; ient < ECS_N_ENTMAIL; ient++) {

    if (cpt_elt_ent[ient] > 0) {

      ECS_MALLOC(elt_pos_som_ent[ient],
                 cpt_elt_ent[ient] + 1,
                 ecs_size_t);

      elt_pos_som_ent[ient][0] = 1;

      ECS_MALLOC(elt_val_som_ent[ient],
                 cpt_val_som_ent[ient],
                 ecs_int_t);

      ECS_MALLOC(elt_val_color_ent[ient]  ,
                 cpt_elt_ent[ient],
                 ecs_int_t);
    }
  }

  /* Remise à zéro des compteurs de dimensionnement */

  for (ient = ECS_ENTMAIL_FAC; ient < ECS_N_ENTMAIL; ient++) {

    cpt_elt_ent[ient]     = 0;
    cpt_val_som_ent[ient] = 0;

  }

  /* Boucle sur la liste des éléments pour construction */
  /* des variables de stockage pour le transfert dans   */
  /* la structure de maillage                           */
  /*----------------------------------------------------*/

  for (ielts_loc = 0; ielts_loc < taille_liste_elems; ielts_loc++) {

    elems_loc = (*liste_elems) + ielts_loc;

    /* Identification de l'entité concernée */
    /*--------------------------------------*/

    typ_geo = elems_loc->elt_typ;
    entmail_e = ecs_maillage_pre__ret_typ_geo(elems_loc->elt_typ);


    /* Couleur des éléments lus positionnée au numéro du part auquel
       appartiennent les éléments */

    for (icoul[entmail_e] = 0;
         icoul[entmail_e] < cpt_coul_ent[entmail_e]
           && val_coul_ent[entmail_e][icoul[entmail_e]] != elems_loc->num_part;
         icoul[entmail_e]++);

    if (icoul[entmail_e] == cpt_coul_ent[entmail_e]) {

      /* La valeur de la couleur n'a pas encore été stockée */

      ECS_REALLOC(val_coul_ent[entmail_e],
                  cpt_coul_ent[entmail_e] + 1,
                  ecs_int_t);
      ECS_REALLOC(cpt_elt_coul_ent[entmail_e],
                  cpt_coul_ent[entmail_e] + 1,
                  ecs_size_t);

      cpt_elt_coul_ent[entmail_e][icoul[entmail_e]] = 0;
      val_coul_ent[entmail_e][icoul[entmail_e]] = elems_loc->num_part;
      cpt_coul_ent[entmail_e]++;
    }

    /* Traitement des éléments "classiques" */
    /*--------------------------------------*/

    if (typ_geo != ECS_ELT_TYP_CEL_POLY && typ_geo != ECS_ELT_TYP_FAC_POLY) {

      nbr_som_elt = ecs_fic_elt_typ_liste_c[typ_geo].nbr_som;

      for (ielt = 0; ielt < elems_loc->nbr_ele; ielt++) {

        /* Construction connectivité */

        pos_elt = elt_pos_som_ent[entmail_e][cpt_elt_ent[entmail_e]];

        elt_pos_som_ent[entmail_e][cpt_elt_ent[entmail_e] + 1]
          = pos_elt + nbr_som_elt;

        for (isom = 0; isom < nbr_som_elt; isom++)
          elt_val_som_ent[entmail_e][pos_elt - 1 + isom]
            = elems_loc->connect[ielt*nbr_som_elt + isom];

        /* Affectation couleurs */

        cpt_elt_coul_ent[entmail_e][icoul[entmail_e]]++;
        elt_val_color_ent[entmail_e][cpt_elt_ent[entmail_e]]
          = icoul[entmail_e] + 1;

        /* Mise à jour compteurs */

        cpt_elt_ent[entmail_e] += 1;
        cpt_val_som_ent[entmail_e] += nbr_som_elt;

      }

    }

    /* Traitement des Polygones */
    /*--------------------------*/

    else if (typ_geo == ECS_ELT_TYP_FAC_POLY) {

      ecs_int_t  cpt_som_loc = 0;

      for (ielt = 0; ielt < elems_loc->nbr_ele; ielt++) {

        nbr_som_elt = elems_loc->nbr_n[ielt];

        /* Construction connectivité */

        pos_elt = elt_pos_som_ent[entmail_e][cpt_elt_ent[entmail_e]];

        elt_pos_som_ent[entmail_e][cpt_elt_ent[entmail_e] + 1]
          = pos_elt + nbr_som_elt;

        for (isom = 0; isom < nbr_som_elt; isom++)
          elt_val_som_ent[entmail_e][pos_elt - 1 + isom]
            = elems_loc->connect[cpt_som_loc++];

        /* Affectation couleurs */

        cpt_elt_coul_ent[entmail_e][icoul[entmail_e]]++;
        elt_val_color_ent[entmail_e][cpt_elt_ent[entmail_e]]
          = icoul[entmail_e] + 1;

        /* Mise à jour compteurs */

        cpt_elt_ent[entmail_e] += 1;
        cpt_val_som_ent[entmail_e] += nbr_som_elt;

      }
    }

    /* Traitement des Polyèdres */
    /*--------------------------*/

    else if (typ_geo == ECS_ELT_TYP_CEL_POLY) {

      ecs_int_t nbr_fac_elt;       /* nombre de faces par élément */
      ecs_int_t ifac_elt;

      ecs_int_t  num_som_deb_fac;
      ecs_int_t  cpt_fac_loc = 0;
      ecs_int_t  cpt_som_loc = 0;

      /* Remplissage de la connectivité */

      for (ielt = 0; ielt < elems_loc->nbr_ele; ielt++) {

        nbr_fac_elt = elems_loc->nbr_f[ielt];

        /* Construction connectivité */

        pos_elt = elt_pos_som_ent[entmail_e][cpt_elt_ent[entmail_e]];

        for (ifac_elt = 0; ifac_elt < nbr_fac_elt; ifac_elt++) {

          nbr_som_fac = elems_loc->nbr_n[cpt_fac_loc];

          /* Le premier noeud de chaque face est répété en queue pour la
             détection de fin de face */

          num_som_deb_fac = elems_loc->connect[cpt_som_loc];

          for (isom = 0; isom < nbr_som_fac; isom++)
            elt_val_som_ent[entmail_e][cpt_val_som_ent[entmail_e] + isom]
              = elems_loc->connect[cpt_som_loc++];

          elt_val_som_ent[entmail_e][cpt_val_som_ent[entmail_e] + isom]
              = num_som_deb_fac;

          /* Mise à jour partielle compteurs */

          cpt_val_som_ent[entmail_e] += nbr_som_fac + 1;

          cpt_fac_loc++;

        }

        elt_pos_som_ent[entmail_e][cpt_elt_ent[entmail_e] + 1]
          = cpt_val_som_ent[entmail_e] + 1;

        /* Affectation couleurs */

        cpt_elt_coul_ent[entmail_e][icoul[entmail_e]]++;
        elt_val_color_ent[entmail_e][cpt_elt_ent[entmail_e]]
          = icoul[entmail_e] + 1;

        /* Mise à jour partielle compteurs */

        cpt_elt_ent[entmail_e] += 1;

      }

    } /* Fin du traitement selon le type d'entité */

    /* Libération mémoire */

    if (elems_loc->nbr_n != NULL)
      ECS_FREE(elems_loc->nbr_n);

    if (elems_loc->nbr_f != NULL)
      ECS_FREE(elems_loc->nbr_f);

    ECS_FREE(elems_loc->connect);

  }

  ECS_FREE(*liste_elems);


  /* Transformation des références en indices */
  /*==========================================*/

  if (liste_noeuds->id_noeud != NULL) {

    for (ient = ECS_ENTMAIL_FAC ; ient < ECS_N_ENTMAIL ; ient++)
      if (elt_val_som_ent[ient] != NULL)
        ecs_maillage_pre__label_en_indice
          (liste_noeuds->nbr_noeuds,
           elt_pos_som_ent[ient][cpt_elt_ent[ient]],
           liste_noeuds->id_noeud,
           elt_val_som_ent[ient]) ;

    ECS_FREE(liste_noeuds->id_noeud) ;
  }

  /* Transfert des valeurs lues dans les structures d'entité de maillage */
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
 *  Détection big-endian/Little-endian au format EnSight 6 binaire C ;
 *  On se place au moment de l'appel au début de la définition des noeuds.
 *----------------------------------------------------------------------------*/

static void
ecs_loc_pre_ens__c_bin_endian_6(ecs_file_t  *fic,
                                bool         lire_id_noeud)
{
  char          chaine[ECS_LOC_LNG_MAX_CHAINE_ENS];

  ecs_int_t     nbr_noeuds;
  int32_t       val_int_lue;

  ecs_int_t     ret = 0;

  bool          swap_endian = false;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (ecs_file_get_type(fic) != ECS_FILE_TYPE_BINARY)
    return;

  ecs_loc_pre_ens__lit_chaine(fic, chaine, NULL);

  if (strncmp(chaine, "coordinates", strlen("coordinates")) != 0)
    ecs_error(__FILE__, __LINE__, 0,
              _("EnSight: field \"%s\" encountered where\n"
                "\"coordinates\" was expected."), chaine);

  ecs_loc_pre_ens__lit_int(fic, 1, &val_int_lue, 8, NULL);

  nbr_noeuds = (ecs_int_t)val_int_lue;

  /* Première permutation big-endian/little-endian si valeur impossible */

  if (nbr_noeuds < 0) {

    swap_endian = true;

    /* On revient en arrière pour la suite */

    ret = ecs_file_seek(fic, -1 * sizeof(int32_t), ECS_FILE_SEEK_CUR);

    if (ret == 0)
      ret = ecs_file_seek(fic, -80 * sizeof(char), ECS_FILE_SEEK_CUR);
  }

  /* Sinon, si la valeur lue pour le nombre de noeuds est plausible,
     on vérifie si elle est effectivement possible */

  else {

    if (lire_id_noeud == true)
      ret = ecs_file_seek(fic,
                          nbr_noeuds * sizeof(int32_t),
                          ECS_FILE_SEEK_CUR);

    if (ret == 0)
      ret = ecs_file_seek(fic,
                          nbr_noeuds * 3 * sizeof(float),
                          ECS_FILE_SEEK_CUR);

    /* Si on arrive effectivement à la fin de fichier, c'est
       que le fichier contient des noeuds, mais pas de "parts" ;
       sinon, on doit avoir un "part", ou c'est que l'on n'a
       pas les bonnes dimensions de tableaux, et donc probablement
       pas la bonne interprétation du nombre de noeuds lus
       -> on doit inverser le paramétrage big-endian/little-endian */

    if (ret == 0 && ecs_file_eof(fic) == 0) {

      if (ecs_file_read_try(chaine, sizeof(char), 80, fic) == 80) {

        if (strncmp(chaine, "part", strlen("part")) != 0)
          swap_endian = true;
      }
      else
        swap_endian = true;
    }

    /* Si l'on n'a pu se positionner, c'est probablement que l'on n'a
       probablement pas la bonne interprétation du nombre de noeuds lus,
       et que l'on a cherché à dépasser la fin du fichier
       -> on doit inverser le paramétrage big-endian/little-endian */

    else if (ret != 0)
      swap_endian = true;

    /* On revient à la position de départ */

    ecs_file_rewind(fic);

    ret = ecs_file_seek(fic, 80 * 5 * sizeof(char), ECS_FILE_SEEK_CUR);
  }

  /* Si l'on n'a pas pu se repositionner, on a une erreur */

  if (ret != 0) {

    ecs_file_read_check_error(fic, 0);

    ecs_error(__FILE__, __LINE__, 0,
              _("EnSight: positioning (seek) error\n"
                "in file \"%s\"."), ecs_file_get_name(fic));
  }

  /* Modification bi-endian/little-endian si nécessaire */

  if (swap_endian == true) {

    if (ecs_file_get_swap_endian(fic) == 1)
      ecs_file_set_swap_endian(fic, 0);
    else
      ecs_file_set_swap_endian(fic, 1);

  }
}

/*----------------------------------------------------------------------------
 *  Lecture d'un fichier au format EnSight Gold
 *   et affectation des données dans la structure de maillage
 *----------------------------------------------------------------------------*/

static ecs_maillage_t  *
ecs_loc_pre_ens__lit_geo_gold(const char       *nom_fic_geo,
                              int               num_part)
{
  char   chaine[ECS_LOC_LNG_MAX_CHAINE_ENS];
  char   chaine_aux[ECS_LOC_LNG_MAX_CHAINE_ENS];

  char    *ret = NULL;

  ecs_int_t   nbr_elt_lus = 0;
  bool        pb_rub;
  bool        lire_node_id;
  bool        lire_elem_id;
  bool        importer_part;

  ecs_file_t  *fic;

  float    coo_min_max[6];

  ecs_int_t        taille_liste_noeuds = 0;
  ecs_int_t        taille_liste_elems = 0;

  ecs_loc_noeuds_ens_t  *liste_noeuds = NULL;
  ecs_loc_elems_ens_t   *liste_elems = NULL;

  int         num_ligne = 1;
  bool        fin_lecture = false;
  bool        affiche_extents = false;

  /* Création d'un maillage initialement vide (valeur de retour) */

  ecs_maillage_t  *maillage = ecs_maillage__cree_nodal();

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Ouverture du fichier geo */
  /*--------------------------*/

  fic = ecs_loc_pre_ens__ouverture_fic_geo(nom_fic_geo);

  /* Lecture de l'entête */
  /*---------------------*/

  /* 2 lignes de description */

  ecs_loc_pre_ens__lit_chaine(fic, chaine, &num_ligne);

  printf(_("\n"
           "  Description:\n\n  %s\n"), chaine);

  ecs_loc_pre_ens__lit_chaine(fic, chaine, &num_ligne);

  printf("  %s\n\n", chaine);

  /* Infos sur ids sommets */

  ecs_loc_pre_ens__lit_chaine(fic, chaine, &num_ligne);

  pb_rub = false;

  if (strncmp(chaine, "node id", strlen("node id")) != 0)
    pb_rub = true;
  else if (sscanf(chaine, "%*s %*s %s", chaine_aux) != 1)
    pb_rub = true;

  if (pb_rub == true)
    ecs_error(__FILE__, __LINE__, 0,
              _("EnSight: \"node id\" field missing or badly placed."));
  else
    printf("  node id :    %s\n", chaine_aux);

  if (   (strncmp(chaine_aux, "given",  strlen("given"))   == 0)
      || (strncmp(chaine_aux, "ignore", strlen("ignore"))  == 0))
    lire_node_id = true;
  else
    lire_node_id = false;

  /* Infos sur ids éléments */

  ecs_loc_pre_ens__lit_chaine(fic, chaine, &num_ligne);

  pb_rub = false;

  if (strncmp(chaine, "element id", strlen("element id")) != 0)
    pb_rub = true;
  else if (sscanf(chaine, "%*s %*s %s", chaine_aux) != 1)
    pb_rub = true;

  if (pb_rub == true)
    ecs_error(__FILE__, __LINE__, 0,
              _("EnSight: \"element id\" field missing or badly placed"));
  else
    printf("  element id : %s\n\n", chaine_aux);

  if (   (strncmp(chaine_aux, "given",  strlen("given"))   == 0)
      || (strncmp(chaine_aux, "ignore", strlen("ignore"))  == 0))
    lire_elem_id = true;
  else
    lire_elem_id = false;

  /* "Extents" ou début de "part" */

  ret = ecs_loc_pre_ens__lit_chaine_essai(fic, chaine, &num_ligne);

  if (ret != NULL && strncmp(chaine, "extents", strlen("extents")) == 0) {

    ecs_loc_pre_ens__lit_float(fic, 6, coo_min_max, &num_ligne);

    /* On ne peut pas encore afficher les "extents", car l'on n'a
       pas encoré détecté l'aspect "big-endian/little-endian"
       d'un fichier binaire C, et les octets peuvent donc être permutés */

    affiche_extents = true;

    ret = ecs_loc_pre_ens__lit_chaine_essai(fic, chaine, &num_ligne);

  }
  else if (ret != NULL && strncmp(chaine, "part", strlen("part")) != 0)
    ecs_error(__FILE__, __LINE__, 0,
              _("EnSight: \"%s\" field encountered where\n"
                "\"extents\" or \"part\" was expected."), chaine);

  if (ret == NULL)
    ecs_error(__FILE__, __LINE__, 0,
              _("File \"%s\"\n"
                "defines no \"part\" (i.e. mesh)."),
              ecs_file_get_name(fic));

  /* Début des "parts" */
  /*-------------------*/

  while (fin_lecture == false)  {

    int32_t  num_part_lu;
    int32_t  cpt_part_lu = 0;

    /* Numéro de part */

    ecs_loc_pre_ens__lit_int(fic, 1, &num_part_lu, 10, &num_ligne);

    /* Détection "big-endian/little-endian" pour un fichier binaire C */

    if (cpt_part_lu == 0 && ecs_file_get_type(fic) == ECS_FILE_TYPE_BINARY) {

      if (num_part_lu < 0 || num_part_lu > 1000) {
        ecs_file_swap_endian(&num_part_lu, &num_part_lu, 4, 1);
        if (num_part_lu < 0 || num_part_lu > 1000)
          ecs_error
            (__FILE__, __LINE__, 0,
             _("EnSight: file \"%s\" seems to be\n"
               "of C binary type but the number of the first \"part\"\n"
               "is < 0 ou > 1000 and provokes the failure of the automatic\n"
               "\"big-endian/little-endian\" detection."),
             ecs_file_get_name(fic));
        else {
          if (ecs_file_get_swap_endian(fic) == 1)
            ecs_file_set_swap_endian(fic, 0);
          else
            ecs_file_set_swap_endian(fic, 1);
          /* Correction extents */
          ecs_file_swap_endian(coo_min_max,
                               coo_min_max,
                               sizeof(float),
                               6);
        }
      }

    }

    /* On peut maintenant afficher les "extents" */

    if (affiche_extents == true) {
      printf(_("  xmin = %12.5e; xmax = %12.5e\n"
               "  ymin = %12.5e; ymax = %12.5e\n"
               "  zmin = %12.5e; zmax = %12.5e\n\n"),
             coo_min_max[0], coo_min_max[1], coo_min_max[2],
             coo_min_max[3], coo_min_max[4], coo_min_max[5]);
      affiche_extents = false;
    }

    if (taille_liste_elems > 0 && lire_node_id == false) {
      printf(_("  Remark: no vertex ids given\n"
               "  --> impossible to merge EnSight \"parts\",\n"
               "      so the following \"parts\" are ignored.\n\n"));
      break;
    }
    else if (num_part > 0 && num_part != num_part_lu) {
      importer_part = false;
      printf(_("  part %2d (ignored): "), (int)num_part_lu);
    }
    else {
      importer_part = true;
      printf(_("  part %2d: "), (int)num_part_lu);
    }

    /* Description */

    ecs_loc_pre_ens__lit_chaine(fic, chaine, &num_ligne);

    printf("%s\n\n", chaine);

    /* Lecture des noeuds */
    /*--------------------*/

    ecs_loc_pre_ens__lit_noeuds_gold(fic,
                                     (ecs_int_t)num_part_lu,
                                     lire_node_id,
                                     importer_part,
                                     &taille_liste_noeuds,
                                     &liste_noeuds,
                                     chaine,
                                     &num_ligne);

    /* Lecture des éléments */
    /*----------------------*/

    do {

      ret = ecs_loc_pre_ens__lit_chaine_essai(fic, chaine, &num_ligne);

      if (ret == NULL)
        fin_lecture = true;

      else

        nbr_elt_lus = ecs_loc_pre_ens__lit_elem_gold
                        (fic,
                         (ecs_int_t)num_part_lu,
                         lire_elem_id,
                         importer_part,
                         liste_noeuds + taille_liste_noeuds - 1,
                         &taille_liste_elems,
                         &liste_elems,
                         chaine,
                         &num_ligne);

    } while (nbr_elt_lus > -1 && fin_lecture == false);

    printf("\n");

    /* Fusion des définitions des noeuds */

    if (importer_part)

      ecs_loc_pre_ens__trie_noeuds(liste_noeuds + taille_liste_noeuds - 1);

    if (taille_liste_noeuds == 2) {

      ecs_loc_pre_ens__fusion_noeuds(liste_noeuds,
                                     liste_noeuds + 1);

      taille_liste_noeuds = 1; /* La réallocation éventuelle à 2 entités
                                   de liste_noeuds[] sera inutile mais
                                   triviale et sans risque */
    }

    /* Autres parts ? */

    cpt_part_lu++;

    if (num_part > 0 && num_part == num_part_lu)
      fin_lecture = true;
  }

  /* Si l'on n'a pas trouvé de "part" (ou pas celle demandée),
     la liste des noeuds n'est encore pas définie */

  if (liste_noeuds == NULL) {
    if (num_part > 0)
      ecs_error(__FILE__, __LINE__, 0,
                _("EnSight: file \"%s\" does not contain a\n"
                  "\"part\" with the required number (%d)."),
                ecs_file_get_name(fic), (int)num_part);
    else
      ecs_error(__FILE__, __LINE__, 0,
                _("EnSight: file \"%s\" does not contain any \"part\"."),
                ecs_file_get_name(fic));
  }

  /* Fermeture du fichier de lecture du maillage */

  ecs_file_free(fic);

  /* Remplissage de la structure de maillage */
  /*-----------------------------------------*/

  /* Traitement des sommets */

  ecs_maillage_pre__cree_som(maillage,
                             liste_noeuds[0].nbr_noeuds,
                             liste_noeuds[0].coord);

  /* Traitement des éléments */

  /* Concaténation des éléments par leur dimension */

  ecs_loc_pre_ens__concat_elems(maillage,
                                taille_liste_elems,
                                &(liste_noeuds[0]),
                                &liste_elems);

  ECS_FREE(liste_noeuds);

  /* Renvoi de la structure de maillage */

  return maillage;
}

/*----------------------------------------------------------------------------
 *  Lecture d'un fichier au format EnSight 6
 *   et affectation des données dans la structure de maillage
 *----------------------------------------------------------------------------*/

static ecs_maillage_t  *
ecs_loc_pre_ens__lit_geo_6(const char       *nom_fic_geo,
                           const ecs_int_t   num_part)
{
  char              chaine[ECS_LOC_LNG_MAX_CHAINE_ENS];
  char              chaine_aux[ECS_LOC_LNG_MAX_CHAINE_ENS];

  char             *ret = NULL;

  ecs_int_t         nbr_elt_lus = 0;
  bool              pb_rub;
  bool              lire_node_id;
  bool              lire_elem_id;
  bool              importer_part;

  ecs_file_t       *fic;

  ecs_int_t         taille_liste_noeuds = 0;
  ecs_int_t         taille_liste_elems = 0;

  ecs_loc_noeuds_ens_t  *liste_noeuds = NULL;
  ecs_loc_elems_ens_t   *liste_elems = NULL;

  int         num_ligne   = 1;
  bool        fin_lecture = false;

  /* Création d'un maillage initialament vide (valeur de retour) */

  ecs_maillage_t  *maillage = ecs_maillage__cree_nodal();

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Ouverture du fichier Géométrie et vérification du format */
  /*----------------------------------------------------------*/

  fic = ecs_loc_pre_ens__ouverture_fic_geo(nom_fic_geo);

  /* Lecture de l'entête */
  /*---------------------*/

  /* 2 lignes de description */

  ecs_loc_pre_ens__lit_chaine(fic, chaine, &num_ligne);

  printf(_("\n"
           "  Description:\n\n  %s\n"), chaine);

  ecs_loc_pre_ens__lit_chaine(fic, chaine, &num_ligne);

  printf("  %s\n\n", chaine);

  /* Infos sur ids sommets */

  ecs_loc_pre_ens__lit_chaine(fic, chaine, &num_ligne);

  pb_rub = false;

  if (strncmp(chaine, "node id", strlen("node id")) != 0)
    pb_rub = true;
  else if (sscanf(chaine, "%*s %*s %s", chaine_aux) != 1)
    pb_rub = true;

  if (pb_rub == true)
    ecs_error(__FILE__, __LINE__, 0,
              _("EnSight: \"node id\" field missing or badly placed."));
  else
    printf("  node id:    %s\n", chaine_aux);

  if (   (strncmp(chaine_aux, "given",  strlen("given"))   == 0)
      || (strncmp(chaine_aux, "ignore", strlen("ignore"))  == 0))
    lire_node_id = true;
  else
    lire_node_id = false;


  /* Infos sur ids éléments */

  ecs_loc_pre_ens__lit_chaine(fic, chaine, &num_ligne);

  pb_rub = false;

  if (strncmp(chaine, "element id", strlen("element id")) != 0)
    pb_rub = true;
  else if (sscanf(chaine, "%*s %*s %s", chaine_aux) != 1)
    pb_rub = true;

  if (pb_rub == true)
    ecs_error(__FILE__, __LINE__, 0,
              _("EnSight: \"element id\" field missing or badly placed"));
  else
    printf("  element id: %s\n\n", chaine_aux);

  if (   (strncmp(chaine_aux, "given",  strlen("given"))   == 0)
      || (strncmp(chaine_aux, "ignore", strlen("ignore"))  == 0))
    lire_elem_id = true;
  else
    lire_elem_id = false;

  /* Détection big-endian/little-endian en cas de fichier binaire C */
  /*----------------------------------------------------------------*/

  if (ecs_file_get_type(fic) == ECS_FILE_TYPE_BINARY)

    ecs_loc_pre_ens__c_bin_endian_6(fic, lire_node_id);


  /* Lecture des noeuds */
  /*--------------------*/

  ecs_loc_pre_ens__lit_noeuds_6(fic,
                                lire_node_id,
                                &taille_liste_noeuds,
                                &liste_noeuds,
                                chaine,
                                &num_ligne);

  ret = ecs_loc_pre_ens__lit_chaine_essai(fic, chaine, &num_ligne);

  if (ret == NULL)
    ecs_error(__FILE__, __LINE__, 0,
              _("File \"%s\"\n"
                "defines no \"part\" (i.e. mesh)."),
              ecs_file_get_name(fic));

  else if (ret != NULL && strncmp(chaine, "part", strlen("part")) != 0)
    ecs_error(__FILE__, __LINE__, 0,
              _("EnSight: \"%s\" field encountered where\n"
                "\"part\" was expected."), chaine);


  /* Début des "parts" */
  /*-------------------*/

  while (fin_lecture == false)  {

    int32_t  num_part_lu;
    int32_t  cpt_part_lu = 0;

    /* Numéro de part */

    num_part_lu = atoi(chaine+strlen("part"));

    /* Détection "big-endian/little-endian" pour un fichier binaire C */


    /* Pour EnSight6 on peut toujours fusionner */

    if (num_part > 0 && num_part != num_part_lu) {
      importer_part = false;
      printf(_("  part %2d (ignored): "), (int)num_part_lu);
    }
    else {
      importer_part = true;
      printf(_("  part %2d: "), (int)num_part_lu);
    }

    /* Description */

    ecs_loc_pre_ens__lit_chaine(fic, chaine, &num_ligne);

    printf("%s\n\n", chaine);

    /* Tri des noeuds */

    ecs_loc_pre_ens__trie_noeuds(liste_noeuds + taille_liste_noeuds - 1);

    /* Lecture des éléments */
    /*----------------------*/

    do {

      ret = ecs_loc_pre_ens__lit_chaine_essai(fic, chaine, &num_ligne);

      if (ret == NULL)
        fin_lecture = true;

      else

        nbr_elt_lus
          = ecs_loc_pre_ens__lit_elem_6(fic,
                                        (ecs_int_t)num_part_lu,
                                        lire_elem_id,
                                        importer_part,
                                        &taille_liste_elems,
                                        &liste_elems,
                                        chaine,
                                        &num_ligne);

    } while (nbr_elt_lus > -1 && fin_lecture == false);

    printf("\n");

    /* Autres parts ? */

    cpt_part_lu++;

    if (num_part > 0 && num_part == num_part_lu)
      fin_lecture = true;

  }

  /* Si l'on n'a pas trouvé de "part" (ou pas celle demandée),
     la liste des noeuds n'est encore pas définie */

  if (liste_elems == NULL) {
    if (num_part > 0)
      ecs_error(__FILE__, __LINE__, 0,
                _("EnSight: file \"%s\" does not contain a\n"
                  "\"part\" with the required number (%d)."),
                ecs_file_get_name(fic), (int)num_part);
    else
      ecs_error(__FILE__, __LINE__, 0,
                _("EnSight: file \"%s\" does not contain any \"part\"."),
                ecs_file_get_name(fic));
  }

  /* Fermeture du fichier de lecture du maillage */

  ecs_file_free(fic);

  /* Remplissage de la structure de maillage */
  /*-----------------------------------------*/

  /* Traitement des sommets */

  ecs_maillage_pre__cree_som(maillage,
                             liste_noeuds[0].nbr_noeuds,
                             liste_noeuds[0].coord);

  /* Traitement des éléments */

  /* Concaténation des éléments par leur dimension */

  ecs_loc_pre_ens__concat_elems(maillage,
                                taille_liste_elems,
                                &(liste_noeuds[0]),
                                &liste_elems);

  ECS_FREE(liste_noeuds);

  /* Renvoi de la structure de maillage */

  return maillage;
}

/*============================================================================
 *                             Fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Lecture d'un fichier au format EnSight Gold
 *   et affectation des données dans la structure de maillage
 *----------------------------------------------------------------------------*/

ecs_maillage_t *
ecs_pre_ens__lit_maillage(const char  *nom_fic_case,
                          int          num_maillage)
{
  char             *nom_fic_geo;
  int               timeset;
  int               fileset;
  size_t            ind;
  ecs_file_t       *fic;

  char              ligne[ECS_LOC_LNG_MAX_CHAINE_ENS];
  char              nom_fic_geo_base[ECS_LOC_LNG_MAX_CHAINE_ENS];

  int               num_ligne = 1;

  bool              fmt_ensight      = false;
  bool              fmt_ensight_gold = false;

  char             *ret = NULL;

  ecs_maillage_t   *maillage = NULL;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  printf(_("\n\n"
           "Reading mesh from file in EnSight format\n"
           "----------------------\n"));

  printf(_("  \"case\" file: %s\n"), nom_fic_case);

  /* Ouverture du fichier Case */

  fic = ecs_file_open(nom_fic_case,
                      ECS_FILE_MODE_READ,
                      ECS_FILE_TYPE_TEXT);

  /* Vérification du format */

  do {
    ret = ecs_file_gets_try(ligne, ECS_LOC_LNG_MAX_CHAINE_ENS, fic, &num_ligne);
  } while (ret != NULL && strncmp(ligne, "FORMAT", strlen("FORMAT")) != 0);

  if (ret != NULL) {

    do {
      ret = ecs_file_gets_try(ligne, ECS_LOC_LNG_MAX_CHAINE_ENS,
                              fic, &num_ligne);
    } while (ret != NULL && strncmp(ligne, "type:", strlen("type:")) != 0);

  }

  if (ret != NULL) {

    for (ind = strlen("type:");
         ligne[ind] != '\0' && (ligne[ind] == ' ' || ligne[ind] == '\t');
         ind++);

    if (strncmp(ligne + ind, "ensight", strlen("ensight")) == 0) {

      fmt_ensight = true;

      ind += strlen("ensight");
      while (ligne[ind] != '\0' && (ligne[ind] == ' ' || ligne[ind] == '\t'))
        ind++;

      if (strncmp(ligne + ind, "gold", strlen("gold")) == 0)

        fmt_ensight_gold = true;

    }

  }

  if (fmt_ensight == false)
    ecs_error(__FILE__, __LINE__, 0,
              _("File \"%s\" does not seem to be a valid\n"
                "EnSight 6 or Gold case file."),
              nom_fic_case);

  /* Recherche des infos sur le fichier géométrique */

  do {
    ret = ecs_file_gets_try(ligne, ECS_LOC_LNG_MAX_CHAINE_ENS,
                            fic, &num_ligne);
  } while (ret != NULL && strncmp(ligne, "GEOMETRY", strlen("GEOMETRY")) != 0);

  if (ret != NULL) {

    do {
      ret = ecs_file_gets_try(ligne, ECS_LOC_LNG_MAX_CHAINE_ENS,
                              fic, &num_ligne);
    } while (ret != NULL && strncmp(ligne, "model:", strlen("model:")) != 0);

  }

  if (ret != NULL) {

    /* La rubrique model: contient deux numéros optionnels (numéro de pas de
       temps et de jeux de fichiers), le nom de base du ou des fichiers
       géométriques, et éventuellement l'option "change_coords_only" */

    if (sscanf(ligne, "%*s %d %d %s",
               &timeset, &fileset, nom_fic_geo_base) != 3) {
      if (sscanf(ligne, "%*s %d %s",
                 &timeset, nom_fic_geo_base) != 2) {
        if (sscanf(ligne, "%*s %s",
                   nom_fic_geo_base) != 1)
          ecs_error(__FILE__, __LINE__, 0,
                    _("The \"%s\" case file does not seem to\n"
                      "indicate a geometry file"),
                    nom_fic_case);
      }
    }

    /* On vérifie que le nom ne contienne pas de caractères "*" */

    for (ind = 0; nom_fic_geo_base[ind] != '\0'; ind++)
      if (nom_fic_geo_base[ind] == '*')
        ecs_error(__FILE__, __LINE__, 0,
                  _("The \"%s\" case file seems to indicate the\n"
                    "series of geometric files named:\n"
                    "\"%s\".\n"
                    "A single file must be chosen."),
                  nom_fic_case, nom_fic_geo_base);

  }

  /* On n'a plus besoin du fichier ".case" */

  ecs_file_free(fic);

  /* Maintenant, on extrait le préfixe du nom du fichier ".case" */

  for (ind = strlen(nom_fic_case) - 1;
       ind > 0 && nom_fic_case[ind] != ECS_PATH_SEP;
       ind--);

  if (nom_fic_case[ind] == ECS_PATH_SEP)
    ind++;

  ECS_MALLOC(nom_fic_geo, ind + strlen(nom_fic_geo_base) + 1, char);
  strncpy(nom_fic_geo, nom_fic_case, ind);
  strcpy(nom_fic_geo + ind, nom_fic_geo_base);

  /* On connaît maintenant le nom du fichier géométrique */

  printf(_("  \"geo\"  file: %s\n\n"), nom_fic_geo);

  /* Lecture du fichier géométrique associé */
  /*----------------------------------------*/

  if (fmt_ensight_gold == true)
    maillage = ecs_loc_pre_ens__lit_geo_gold(nom_fic_geo,
                                             num_maillage);

  else if (fmt_ensight == true)
    maillage = ecs_loc_pre_ens__lit_geo_6(nom_fic_geo,
                                          num_maillage);

  /* Libération mémoire et retour */

  ECS_FREE(nom_fic_geo);

  return maillage;
}

/*----------------------------------------------------------------------------*/
