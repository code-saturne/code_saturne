/*============================================================================
 *  Définition de la fonction
 *   de lecture d'un méta-fichier de maillage
 *============================================================================*/

/*
  This file is part of the Code_Saturne Preprocessor, element of the
  Code_Saturne CFD tool.

  Copyright (C) 1999-2009 EDF S.A., France

  contact: saturne-support@edf.fr

  The Code_Saturne Preprocessor is free software; you can redistribute it
  and/or modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  The Code_Saturne Preprocessor is distributed in the hope that it will be
  useful, but WITHOUT ANY WARRANTY; without even the implied warranty
  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with the Code_Saturne Preprocessor; if not, write to the
  Free Software Foundation, Inc.,
  51 Franklin St, Fifth Floor,
  Boston, MA  02110-1301  USA
*/


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


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_pre.h"


/*----------------------------------------------------------------------------
 *  Fichier  `include' du  paquetage courant associe au fichier courant
 *----------------------------------------------------------------------------*/

#include "ecs_pre_meta.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' privés   du  paquetage courant
 *----------------------------------------------------------------------------*/


/*============================================================================
 *                    Déclaration de paramètres et macros
 *============================================================================*/

/* Pour une lecture de 512 caracteres par ligne  */
/* auxquels il faut ajouter le `\n' et le `\0'  */
/* pour l'affectation dans la chaîne réceptrice */
#define ECS_LOC_LNG_MAX_CHAINE_META  514

/*============================================================================
 *                  Définition de structures locales
 *============================================================================*/


/*============================================================================
 *                       Prototypes de fonctions privées
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Lecture d'un fichier selon une ligne de description
 *   et affectation des données dans une structure de maillage
 *----------------------------------------------------------------------------*/

static ecs_maillage_t  *
ecs_loc_pre_meta__lit_maillage(char        *ligne,
                               const char  *nom_fic,
                               int          num_ligne);

/*============================================================================
 *                             Fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Lecture d'un méta-fichier de maillage
 *   et affectation des données dans la structure de maillage
 *----------------------------------------------------------------------------*/

ecs_maillage_t *
ecs_pre_meta__lit_maillage(const char  *nom_fic)
{
  size_t             ind, ind_tmp, len;
  ecs_file_t       * fic = NULL;

  char               ligne[ECS_LOC_LNG_MAX_CHAINE_META];

  bool               par_transf = false;
  double             transf[3][4] = {{1., 0., 0., 0.},
                                     {0., 1., 0., 0.},
                                     {0., 0., 1., 0.}};

  int                num_ligne = 0;

  ecs_maillage_t   * maillage = NULL;
  ecs_maillage_t   * maillage_tmp = NULL;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  printf(_("\n\n"
           "Reading mesh based on meta-mesh file\n"
           "------------------------------------\n"));

  printf(_("  Meta-mesh file: %s\n\n\n"),
         nom_fic);


  /* Ouverture du fichier Case */

  fic = ecs_file_open(nom_fic,
                      ECS_FILE_MODE_READ,
                      ECS_FILE_TYPE_TEXT);


  while (   ecs_file_gets_try(ligne,
                              ECS_LOC_LNG_MAX_CHAINE_META,
                              fic,
                              &num_ligne)
         != NULL) {

    /* Suppression des commentaires */

    len = strlen(ligne);

    for (ind = 0; ind < len; ind++) {
      if (ligne[ind] == '#') {
        ligne[ind] = '\0';
        len = ind;
        break;
      }
    }

    /* Remplacement des tabulations et par des blancs */

    for (ind = 0; ind < len; ind++) {
      if (ligne[ind] == '\t')
        ligne[ind] = ' ';
    }

    /* Suppression des blancs en fin et début de ligne */

    for (ind = len; ind > 0; ind--) {
      if (ligne[ind-1] != ' ' && ligne[ind-1] != '\n' && ligne[ind-1] != '\r')
        break;
    }
    if (ind != len) {
      ligne[ind] = '\0';
      len = ind;
    }

    for (ind = 0; ind < len; ind++) {
      if (ligne[ind] != ' ')
        break;
    }
    if (ind > 0) {
      len -= ind;
      memmove(ligne, ligne+ind, len + 1);
    }

    /* Supprimer les blancs multiples à l'intérieur */

    ind_tmp = 1;
    for (ind = 1; ind < len; ind++) {
      if (ligne[ind] != ' ' || ligne[ind_tmp-1] != ' ')
        ligne[ind_tmp++] = ligne[ind];
    }
    ligne[ind_tmp] = '\0';

    /* Ignorer les lignes vides après nettoyage */

    if (len == 0)
      continue;

    /* Interpréter */

    if (strncmp(ligne, "read_mesh:", strlen("read_mesh:")) == 0) {

      maillage_tmp
        = ecs_loc_pre_meta__lit_maillage(ligne + strlen("read_mesh:"),
                                         nom_fic,
                                         num_ligne);

      if (maillage == NULL)
        maillage = maillage_tmp;

      else
        ecs_maillage__concatene_nodal(maillage, maillage_tmp);

    }

    else if (strncmp(ligne,
                     "transformation_matrix:",
                     strlen("transformation_matrix:")) == 0) {

      int i, j;
      int nbr_par_transf = 0;
      double transf_tmp[12];

      while (nbr_par_transf < 12) {
        double *tp = transf_tmp + nbr_par_transf;
        nbr_par_transf
          += sscanf(ligne + strlen("transformation_matrix:"),
                    "%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg",
                    tp, tp+1, tp+2, tp+3,
                    tp+4, tp+5, tp+6, tp+7,
                    tp+8, tp+9, tp+10, tp+11);
        if (nbr_par_transf < 12)
          ecs_file_gets(ligne,
                        ECS_LOC_LNG_MAX_CHAINE_META,
                        fic,
                        &num_ligne);
      }

      for (i = 0; i < 3; i++) {
        for (j = 0; j < 4; j++)
          transf[i][j] = transf_tmp[i*4+j];
      }

      par_transf = true;

    }

  }

  if (par_transf == true)
    ecs_maillage__transf_coo(maillage, (void *)transf);

  /* On a terminé */

  ecs_file_free(fic);

  return maillage;
}

/*============================================================================
 *                              Fonctions privées
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Lecture d'un fichier selon une ligne de description
 *   et affectation des données dans une structure de maillage
 *----------------------------------------------------------------------------*/

static ecs_maillage_t *
ecs_loc_pre_meta__lit_maillage(char        *ligne,
                               const char  *nom_fic,
                               int          num_ligne)
{
  size_t             ind_deb, ind_fin, len;

  ecs_pre_format_t   format = ECS_PRE_FORMAT_NUL;
  int                num_maillage = 1;
  bool               cree_grp_cel_section = false;
  bool               cree_grp_cel_zone = false;
  bool               cree_grp_fac_section = false;
  bool               cree_grp_fac_zone = false;

  bool               lire_format = false;
  bool               lire_num = false;
  bool               lire_grp_cel = false;
  bool               lire_grp_fac = false;
  bool               syntax_error = false;

  const char        *nom_fic_maillage = NULL;
  ecs_maillage_t    *maillage = NULL;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  ind_deb = 0;

  /* Interprétation du nom du fichier */

  len = strlen(ligne);

  if (ligne[ind_deb] == ' ')
    ind_deb++;

  if (ligne[ind_deb] == '"') {
    for (ind_fin = ind_deb+1;
         ligne[ind_fin] != '"' && ligne[ind_fin] != '\0';
         ind_fin++);
    if (ligne[ind_fin] == '"')
      ligne[ind_fin] = '\0';
    nom_fic_maillage = ligne + ind_deb + 1;
  }
  else {
    for (ind_fin = ind_deb+1;
         ligne[ind_fin] != ' ' && ligne[ind_fin] != '\0';
         ind_fin++);
    if (ligne[ind_fin] == ' ')
      ligne[ind_fin] = '\0';
    nom_fic_maillage = ligne + ind_deb;

  }

  /* Options supplémentaires */

  for (ind_deb = ind_fin + 1;
       ind_deb < len && ligne[ind_deb] == ' ';
       ind_deb++);

  if (ind_deb < len) {

    size_t i1, i2;

    for (i1 = ind_deb; i1 < len; i1++)
      if (ligne[i1] == ',' ||ligne[i1] == ';')
        ligne[i1] = ' ';

    while (ind_deb < len) {

      /* Repérage et analyse des rubriques */

      i1 = ind_deb;
      for (i2 = i1;
           i2 < len && ligne[i2] != '=' && ligne[i2] != ' ';
           i2++);

      /* Nom option */

      if (strncmp(ligne+i1, "format", i2-i1) == 0)
        lire_format = true;

      else if (strncmp(ligne+i1, "num", i2-i1) == 0)
        lire_num = true;

      else if (strncmp(ligne+i1, "grp_cel", i2-i1) == 0)
        lire_grp_cel = true;

      else if (strncmp(ligne+i1, "grp_fac", i2-i1) == 0)
        lire_grp_fac = true;

      /* Valeur option */

      if (   lire_format == true || lire_num == true
          || lire_grp_cel == true || lire_grp_fac == true) {

        for (i1 = i2;
             i1 < len && (ligne[i1] == '=' || ligne[i1] == ' ');
             i1++);
        for (i2 = i1;
             i2 < len && ligne[i2] != ' ';
             i2++);
        ligne[i2] = '\0';

      }

      if (lire_format == true) {
        format = ecs_pre__type_format(nom_fic_maillage,
                                      ligne + i1);
        lire_format = false;
      }

      else if (lire_num == true) {
        if (sscanf(ligne+i1, "%d", &num_maillage) != 1)
          syntax_error = true;
        lire_num = false;
      }

      else if (lire_grp_cel == true) {
        if (strncmp(ligne+i1, "section", i2-i1) == 0)
          cree_grp_cel_section = true;
        else if (strncmp(ligne+i1, "zone", i2-i1) == 0)
          cree_grp_cel_zone = true;
        lire_grp_cel = false;
      }

      else if (lire_grp_fac == true) {
        if (strncmp(ligne+i1, "section", i2-i1) == 0)
          cree_grp_cel_section = true;
        else if (strncmp(ligne+i1, "zone", i2-i1) == 0)
          cree_grp_cel_zone = true;
        lire_grp_fac = false;
      }

      /* Passage à l'option suivante */

      for (ind_deb = i2 + 1;
           ind_deb < len && ligne[ind_deb] == ' ';
           ind_deb++);

    }
  }

  if (syntax_error == true)
    ecs_error(__FILE__, __LINE__, 0,
              _("Error in mesh options specification on line %d "
                "of file \"%s\"\n\n"
                "\"%s\" should use a syntax of the form:\n"
                "opt1 = value1, opt2 = value2, ..."),
              num_ligne, nom_fic);

  if (format == ECS_PRE_FORMAT_NUL)
    format = ecs_pre__type_format(nom_fic_maillage, NULL);


  /* Lecture effective du maillage */

  if (format == ECS_PRE_FORMAT_META)

    maillage = ecs_pre_meta__lit_maillage(nom_fic_maillage);

  else

    maillage = ecs_pre__lit_maillage(nom_fic_maillage,
                                     format,
                                     num_maillage,
                                     cree_grp_cel_section,
                                     cree_grp_cel_zone,
                                     cree_grp_fac_section,
                                     cree_grp_fac_zone);

  return maillage;
}

/*----------------------------------------------------------------------------*/

