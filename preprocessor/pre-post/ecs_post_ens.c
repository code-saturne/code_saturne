/*============================================================================
 *  Définitions des fonctions de base
 *   réalisant les sorties pour post-traitement Ensight
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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
#include <ctype.h>  /* toupper() */
#include <stdio.h>
#include <string.h>


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/

#include "ecs_def.h"
#include "ecs_file.h"
#include "ecs_mem.h"


/*----------------------------------------------------------------------------
 *  Fichier `include' du  paquetage courant associé au fichier courant
 *----------------------------------------------------------------------------*/

#include "ecs_post_ens.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' privés du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_post_ens_priv.h"


/*============================================================================
 *                              Fonctions privées
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction écrivant un fichier Case
 *
 *  Cette fonction crée le répertoire correspondant si nécessaire
 *----------------------------------------------------------------------------*/

static void
ecs_loc_post_ens__ecr_case(ecs_post_ens_t  *cas_ens)
{
  ecs_file_t  *fic_case;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (cas_ens->nom_fic_case == NULL) {

    size_t ind;
    char *nom_base;

    ECS_MALLOC(nom_base,
               strlen(cas_ens->nom_cas) + 1 + strlen(".ensight") + 1,
               char);
    strcpy(nom_base, cas_ens->nom_cas);
    strcat(nom_base, ".ensight");

    for (ind = 0; ind < strlen(cas_ens->nom_cas); ind++) {
      if (nom_base[ind] == ' ')
        nom_base[ind] = '_';
      else
        nom_base[ind] = tolower(nom_base[ind]);
    }

    /* Création du répertoire si possible */
    /*------------------------------------*/

    if (ecs_file_mkdir_default(nom_base) == 0) {

      /* Nom du répertoire */

      ECS_MALLOC(cas_ens->prefixe_rep,
                 strlen(nom_base) + strlen("/") + 1,
                 char);

      strcpy(cas_ens->prefixe_rep, nom_base);
      strcat(cas_ens->prefixe_rep, "/");

        /* Préfixe des noms de fichiers */

        ECS_MALLOC(cas_ens->prefixe_fic, strlen("mesh") + 1, char);

        strcpy(cas_ens->prefixe_fic, "mesh");

        /* Cas du fichier ".case" (nom en majuscules) */

        ECS_MALLOC(cas_ens->nom_fic_case,
                   strlen(cas_ens->prefixe_rep) + strlen(cas_ens->prefixe_fic)
                   + strlen(".case") + 1, char);

        strcpy(cas_ens->nom_fic_case, cas_ens->prefixe_rep);
        strcat(cas_ens->nom_fic_case, cas_ens->prefixe_fic);

        for (ind = strlen(cas_ens->prefixe_rep);
             cas_ens->nom_fic_case[ind] != '\0';
             ind++)
          cas_ens->nom_fic_case[ind] = toupper(cas_ens->nom_fic_case[ind]);

        strcat(cas_ens->nom_fic_case, ".case");

    }
    else {

      /* Nom du répertoire */

      ECS_MALLOC(cas_ens->prefixe_rep, 1, char);

      strcpy(cas_ens->prefixe_rep, "");

      /* Préfixe des noms de fichiers */

      ECS_MALLOC(cas_ens->prefixe_fic,
                 strlen(nom_base) + strlen("_") + 1,
                 char);

      strcpy(cas_ens->prefixe_fic, nom_base);
      strcat(cas_ens->prefixe_fic, "_");

      /* Cas du fichier ".case" (nom en majuscules) */

      ECS_MALLOC(cas_ens->nom_fic_case,
                 strlen(nom_base) + strlen(".case") + 1,
                 char);

      for (ind = 0; nom_base[ind] != '\0'; ind++)
        nom_base[ind] = toupper(nom_base[ind]);

      strcpy(cas_ens->nom_fic_case, nom_base);
      strcat(cas_ens->nom_fic_case, "_");

      strcat(cas_ens->nom_fic_case, ".case");
    }

    ECS_FREE(nom_base);

    /* Info sur la création du fichier Case */
    /*--------------------------------------*/

    printf("  %s %s\n", _("Creating file:"), cas_ens->nom_fic_case);

  }

  /* Ouverture du fichier Case */
  /*---------------------------*/

  fic_case = ecs_file_open(cas_ens->nom_fic_case,
                           ECS_FILE_MODE_WRITE,
                           ECS_FILE_TYPE_TEXT);

  /* Écriture de la section FORMAT */
  /*-------------------------------*/

  ecs_file_printf(fic_case,
                  "FORMAT\n"
                  "type: ensight gold\n\n");

  /* Écriture de la section GEOMETRY */
  /*---------------------------------*/

  ecs_file_printf(fic_case, "GEOMETRY\n");

  ecs_file_printf(fic_case,
                  "model: %s.geo\n\n",
                  cas_ens->prefixe_fic);

  /* Fermeture du fichier Case */
  /*---------------------------*/

  ecs_file_free(fic_case);
  cas_ens->modifie = false;
}

/*----------------------------------------------------------------------------
 *  Fonction vidant une structure `ecs_post_ens_part_t`
 *----------------------------------------------------------------------------*/

static ecs_post_ens_part_t  *
ecs_loc_post_ens__detruit_part(ecs_post_ens_part_t  *this_part)
{
  ecs_int_t ind;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  ECS_FREE(this_part->nom_part);

  this_part->num_part = 0;
  this_part->nbr_som  = 0;

  if (this_part->nbr_typ_ele > 0) {

    for (ind = 0; ind < this_part->nbr_typ_ele; ind ++)
      ECS_FREE(this_part->nom_typ_ele[ind]);

    ECS_FREE(this_part->nbr_ele_typ);
    ECS_FREE(this_part->nom_typ_ele);

  }

  if (this_part->lst_parents != NULL)
    ECS_FREE(this_part->lst_parents);

  ECS_FREE(this_part);

  return NULL;
}

/*============================================================================
 *                             Fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction initialisant une structure `ecs_post_ens_t`.
 *----------------------------------------------------------------------------*/

ecs_post_ens_t  *
ecs_post_ens__cree_cas(const char  *nom_cas)
{
  ecs_post_ens_t  *cas_ens;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Création de la structure `ecs_post_ens_t` */
  /*-------------------------------------------*/

  ECS_MALLOC(cas_ens, 1, ecs_post_ens_t);

  /* Construction du nom */
  /*---------------------*/

  ECS_MALLOC(cas_ens->nom_cas, strlen(nom_cas) + 1, char);
  strcpy(cas_ens->nom_cas, nom_cas);

  /* Noms qui seront initialisés à la première écriture du fichier case */

  cas_ens->prefixe_rep = NULL;
  cas_ens->prefixe_fic = NULL;
  cas_ens->nom_fic_case = NULL;

  /* Champs relatifs au fichier géométrie */
  /*--------------------------------------*/

  cas_ens->nbr_part = 0;
  cas_ens->tab_part = NULL;
  cas_ens->fic_geo  = NULL;

  cas_ens->modifie = false;

  return cas_ens;
}

/*----------------------------------------------------------------------------
 *  Fonction détruisant une structure `ecs_post_ens_t`.
 *----------------------------------------------------------------------------*/

ecs_post_ens_t  *
ecs_post_ens__detruit_cas(ecs_post_ens_t  *cas_ens)
{
  ecs_int_t ind;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (cas_ens != NULL) {

    /* Destruction du cas */

    ECS_FREE(cas_ens->nom_cas);

    ECS_FREE(cas_ens->prefixe_rep);
    ECS_FREE(cas_ens->prefixe_fic);
    ECS_FREE(cas_ens->nom_fic_case);

    for (ind = 0; ind < cas_ens->nbr_part; ind ++)
      ecs_loc_post_ens__detruit_part(cas_ens->tab_part[ind]);

    if (cas_ens->nbr_part != 0)
      ECS_FREE(cas_ens->tab_part);

    if (cas_ens->fic_geo != NULL)
      ecs_file_free(cas_ens->fic_geo);

    ECS_FREE(cas_ens);
  }

  return cas_ens;
}

/*----------------------------------------------------------------------------
 *  Écriture d'une chaîne de caractères dans un fichier EnSight Gold.
 *----------------------------------------------------------------------------*/

void
ecs_post_ens__ecr_chaine(const ecs_file_t  *fic,
                         const char        *chaine)
{
  char ligne[83];
  size_t ind;

  strncpy(ligne, chaine, 80);

  for (ind = strlen(chaine); ind < 80; ind++)
    ligne[ind] = ' ';
  ecs_file_write(ligne, sizeof(char), 80, fic);
}

/*----------------------------------------------------------------------------
 *  Écriture d'un entier dans un fichier EnSight Gold.
 *----------------------------------------------------------------------------*/

void
ecs_post_ens__ecr_int(const ecs_file_t  *fic,
                      int                val)
{
  int32_t _val = val;

  ecs_file_write(&_val, sizeof(int32_t), 1, fic);
}

/*----------------------------------------------------------------------------
 *  Fonction écrivant le fichier contenant la géométrie
 *----------------------------------------------------------------------------*/

ecs_file_t  *
ecs_post_ens__ecrit_fic_geo(ecs_post_ens_t  *cas_ens)
{
  char   ligne_cas[81];
  char  *nom_cas;
  char  *nom_fic_geo;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Si le répertoire et le fichier .case n'ont pas encore été
     crées, on les crée. */

  ecs_loc_post_ens__ecr_case(cas_ens);

  /* Si le fichier géométrie est simplement fermé, on le réouvre */
  /*-------------------------------------------------------------*/

  if (cas_ens->fic_geo != NULL) {
    ecs_file_open_stream(cas_ens->fic_geo, ECS_FILE_MODE_APPEND);
    return cas_ens->fic_geo;
  }

  /* Construction du nom du fichier contenant la géométrie */
  /*-------------------------------------------------------*/

  ECS_MALLOC(nom_fic_geo,
             strlen(cas_ens->prefixe_rep) + strlen(cas_ens->prefixe_fic)
             + strlen(".geo") + 1,
             char);

  strcpy(nom_fic_geo, cas_ens->prefixe_rep);
  strcat(nom_fic_geo, cas_ens->prefixe_fic);
  strcat(nom_fic_geo, ".geo");

  /* Ouverture du fichier contenant la géométrie */
  /*---------------------------------------------*/

  cas_ens->fic_geo = ecs_file_open(nom_fic_geo,
                                   ECS_FILE_MODE_WRITE,
                                   ECS_FILE_TYPE_BINARY);

  ecs_post_ens__ecr_chaine(cas_ens->fic_geo, "C binary");

  /* Info sur la création du fichier contenant la géométrie */
  /*--------------------------------------------------------*/

  printf("  %s %s\n", _("Creating file:"), nom_fic_geo);

  /* Écriture des 2 premières lignes de commentaires */
  /*-------------------------------------------------*/

  ECS_MALLOC(nom_cas, strlen(cas_ens->prefixe_rep) + 1, char);
  strcpy(nom_cas, cas_ens->prefixe_rep);
  strtok(nom_cas, ".");

  ecs_post_ens__ecr_chaine(cas_ens->fic_geo,
                           "EnSight Gold output by Code_Saturne Preprocessor");

  strcpy(ligne_cas, "Case name: ");
  strncpy(ligne_cas + strlen(ligne_cas), nom_cas, 80 - strlen(ligne_cas));

  ecs_post_ens__ecr_chaine(cas_ens->fic_geo, ligne_cas);

  ECS_FREE(nom_cas);

  ECS_FREE(nom_fic_geo);

  /* Ensight se chargera d'affecter des labels */
  /*-------------------------------------------*/

  ecs_post_ens__ecr_chaine(cas_ens->fic_geo, "node id assign");
  ecs_post_ens__ecr_chaine(cas_ens->fic_geo, "element id assign");

  return cas_ens->fic_geo;
}

/*----------------------------------------------------------------------------*/


