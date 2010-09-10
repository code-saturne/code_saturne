/*============================================================================
 *  Définitions des fonctions de base
 *   réalisant les sorties pour post-traitement Ensight
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
 *  Fonction qui met à jour la structure ecs_post_ens_t en fonction d'une
 *  nouvelle variable : si le temps physique ou la variable correspondante
 *  n'est pas présent dans la structure, on ajoute les elements
 *  correspondants à la structure.
 *----------------------------------------------------------------------------*/

static char *
ecs_loc_post_ens__ajout_var(ecs_post_ens_t  *cas_ens,
                            const char      *nom_champ,
                            ecs_int_t       *num_sor_prec)
{
  size_t      ind;
  char       *nom_var;
  size_t      lng_nom_var;
  size_t      lng_ligne_var;
  ecs_int_t   ivar;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* On écrit la ligne contenant la référence de la variable                  */
  /*  dans le fichier `.case' sous la forme suivante :                        */
  /*                                                                          */
  /* <nom_rubrique> <nom_var> <prefixe>.<nom_var>                             */
  /*                                                                          */
  /* Or la ligne ne doit pas dépasser 79 caractères                           */
  /*  d'après les spécifications d'Ensight.                                   */
  /* Le nom de la variable écrit sur fichier Case d'Ensight est donc          */
  /*  le nom du champ, raccourcit éventuellement pour ne dépasser             */
  /*  les 79 caractères autorises sur la ligne                                */
  /*                                                                          */
  /* `scalar per element:' (19 caractères)                                    */

  lng_ligne_var = 19 + 1 + 1 + strlen(cas_ens->prefixe_fic) + 1;

  lng_nom_var = 79 - lng_ligne_var;

  /* Le nom de la variable est aussi utilisée comme suffixe du fichier : */

  lng_nom_var /= 2; /* Division entière */

  lng_nom_var = ECS_MIN(lng_nom_var, strlen(nom_champ));

  ECS_MALLOC(nom_var, lng_nom_var + 1, char);

  sprintf(nom_var, "%*.*s", (int)lng_nom_var, (int)lng_nom_var, nom_champ);

  for (ind = 0; ind < lng_nom_var; ind++) {
    if (nom_var[ind] == '-')
      nom_var[ind] = '_';
  }

  /* On cherche si cette variable existe déjà  */

  ivar = 0;

  while (ivar < cas_ens->nbr_var &&
         strcmp(cas_ens->nom_var[ivar], nom_var) != 0)
      ivar ++;

  /* Si la variable n'existe pas encore, ivar = cas_ens->nbr_var   */

  if (ivar == cas_ens->nbr_var) {

    /* ----------------------------------------------------------- */
    /* Mise à jour de la structure si la ligne n'existe pas encore */
    /* ----------------------------------------------------------- */

    size_t    lng_imp_nom_var;

    /* Mise à jour du nombre de variables (on se base plutôt sur */
    /* la variable locale ivar à l'intérieur de cette fonction)  */

    cas_ens->nbr_var += 1;

    /* Ajout de nom_var */

    ECS_REALLOC(cas_ens->nom_var, ivar + 1, char *);

    ECS_MALLOC(cas_ens->nom_var[ivar], strlen(nom_var) + 1, char);

    strcpy(cas_ens->nom_var[ivar], nom_var);

    /* Ajout de ligne_var (on choisit de fournir au moins 24 caractères
       pour la zone d'écriture du nom de la variable) */

    ECS_REALLOC(cas_ens->ligne_var, ivar + 1, char *);

    lng_imp_nom_var = ECS_MAX(lng_nom_var, 24);

    ECS_MALLOC(cas_ens->ligne_var[ivar],
               lng_ligne_var + lng_imp_nom_var + lng_nom_var + 1,
               char);

    /* Variable pas encore définie */

    *num_sor_prec = -1;

    /* Nom dans le fichier case */

    sprintf(*(cas_ens->ligne_var + ivar),
            "scalar per element: %-*s %s.%s",
            (int)lng_imp_nom_var,
            nom_var,
            cas_ens->prefixe_fic,
            nom_var);

    cas_ens->modifie = true;

  }
  else {

    *num_sor_prec = 0;

  }

  /* On renvoie le nom associé à la variable */

  return nom_var;
}

/*----------------------------------------------------------------------------
 *  Fonction écrivant un fichier Case
 *----------------------------------------------------------------------------*/

static void
ecs_loc_post_ens__ecr_case(ecs_post_ens_t  *cas_ens)
{
  ecs_file_t    * fic_case;
  ecs_int_t       ivar;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

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

  /* Écriture de la section VARIABLE */
  /*---------------------------------*/

  if (cas_ens->nbr_var > 0) {
    ecs_file_printf(fic_case, "VARIABLE\n");

    for (ivar = 0; ivar < cas_ens->nbr_var; ivar ++) {

      ecs_file_printf(fic_case, "%s\n", cas_ens->ligne_var[ivar]);

    }
  }

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
ecs_post_ens__cree_cas(const char  *nom_cas,
                       bool         no_poly,
                       bool         text,
                       bool         big_endian)
{
  char            *nom_base;
  ecs_post_ens_t  *cas_ens;

  size_t           ind;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Création de la structure `ecs_post_ens_t` */
  /*-------------------------------------------*/

  ECS_MALLOC(cas_ens, 1, ecs_post_ens_t);

  cas_ens->no_poly = no_poly;

  /* Construction du nom */
  /*---------------------*/

  ECS_MALLOC(cas_ens->nom_cas, strlen(nom_cas) + 1, char);
  strcpy(cas_ens->nom_cas, nom_cas);

  ECS_MALLOC(nom_base, strlen(nom_cas) + 1 + strlen(".ensight") + 1, char);
  strcpy(nom_base, nom_cas);
  strcat(nom_base, ".ensight");

  for (ind = 0; ind < strlen(nom_cas); ind++) {
    if (nom_base[ind] == ' ')
      nom_base[ind] = '_';
    else
      nom_base[ind] = tolower(nom_base[ind]);
  }

  /* Création du répertoire si possible */
  /*------------------------------------*/

  if (ecs_file_mkdir_default(nom_base) == 0) {

    /* Nom du répertoire */

    ECS_MALLOC(cas_ens->prefixe_rep, strlen(nom_base) + strlen("/") + 1, char);

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

    ECS_MALLOC(cas_ens->prefixe_fic, strlen(nom_base) + strlen("_") + 1, char);

    strcpy(cas_ens->prefixe_fic, nom_base);
    strcat(cas_ens->prefixe_fic, "_");

    /* Cas du fichier ".case" (nom en majuscules) */

    ECS_MALLOC(cas_ens->nom_fic_case, strlen(nom_base)
               + strlen(".case") + 1, char);

    for (ind = 0; nom_base[ind] != '\0'; ind++)
      nom_base[ind] = toupper(nom_base[ind]);

    strcpy(cas_ens->nom_fic_case, nom_base);
    strcat(cas_ens->nom_fic_case, "_");

    strcat(cas_ens->nom_fic_case, ".case");
  }

  ECS_FREE(nom_base);

  /* Autres membres de la structure `ecs_post_ens_t` */
  /*-------------------------------------------------*/

  cas_ens->nbr_var      = 0;

  cas_ens->nom_var     = NULL;
  cas_ens->ligne_var   = NULL;

  /* Champs relatifs au fichier géométrie */
  /*--------------------------------------*/

  cas_ens->nbr_part = 0;
  cas_ens->tab_part = NULL;
  cas_ens->fic_geo  = NULL;

  /* Première écriture du fichier Case */
  /*-----------------------------------*/

  cas_ens->modifie = true;

  ecs_loc_post_ens__ecr_case(cas_ens);

  cas_ens->modifie = false;

  cas_ens->text = text;
  cas_ens->big_endian = big_endian;

  /* Info sur la création du fichier Case */
  /*--------------------------------------*/

  printf("  %s %s\n", _("Creating file:"), cas_ens->nom_fic_case);

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

    if (cas_ens->nbr_var > 0) {

      for (ind = 0; ind < cas_ens->nbr_var; ind ++)
        ECS_FREE(cas_ens->nom_var[ind]);

      ECS_FREE(cas_ens->nom_var);

      for (ind = 0; ind < cas_ens->nbr_var; ind ++)
        ECS_FREE(cas_ens->ligne_var[ind]);

      ECS_FREE(cas_ens->ligne_var);
    }

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

  if (ecs_file_get_type(fic) != ECS_FILE_TYPE_TEXT) {

    for (ind = strlen(chaine); ind < 80; ind++)
      ligne[ind] = ' ';
    ecs_file_write(ligne, sizeof(char), 80, fic);

  }
  else {

    ligne[80] = '\0';
    ecs_file_printf(fic, "%s\n", ligne);

  }
}

/*----------------------------------------------------------------------------
 *  Écriture d'un entier dans un fichier EnSight Gold.
 *----------------------------------------------------------------------------*/

void
ecs_post_ens__ecr_int(const ecs_file_t  *fic,
                      int                val)
{
  ecs_int_32_t _val;

  if (ecs_file_get_type(fic) != ECS_FILE_TYPE_TEXT) {

    _val = val;
    ecs_file_write(&_val, sizeof(ecs_int_32_t), 1, fic);

  }
  else {

    ecs_file_printf(fic, "%10d\n", val);

  }
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

  if (cas_ens->text == false) {

    cas_ens->fic_geo = ecs_file_open(nom_fic_geo,
                                     ECS_FILE_MODE_WRITE,
                                     ECS_FILE_TYPE_BINARY);

    if (cas_ens->big_endian == true)
      ecs_file_set_big_endian(cas_ens->fic_geo);

    ecs_post_ens__ecr_chaine(cas_ens->fic_geo, "C binary");

  }
  else
    cas_ens->fic_geo = ecs_file_open(nom_fic_geo,
                                     ECS_FILE_MODE_WRITE,
                                     ECS_FILE_TYPE_TEXT);

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

/*----------------------------------------------------------------------------
 *  Fonction construisant le descripteur (`ecs_file_t') du fichier
 *  contenant les valeurs de la variable   sortir pour post-traitement Ensight
 *
 *  La fonction détermine aussi la ligne spécifiant la variable
 *   devant figurer dans le fichier Case
 *----------------------------------------------------------------------------*/

ecs_file_t  *
ecs_post_ens__ecrit_fic_var(ecs_post_ens_t  *cas_ens,
                            const char      *nom_champ)
{
  ecs_file_t       * fic_var;
  ecs_file_mode_t    mode_ouverture;
  char             * nom_var;
  char             * nom_fic_var;
  char               buf[81];

  ecs_int_t          num_sor_prec;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Construction du nom de la variable à partir du nom du champ */
  /*-------------------------------------------------------------*/

  nom_var = ecs_loc_post_ens__ajout_var(cas_ens,
                                        nom_champ,
                                        &num_sor_prec);

  /* Création du descripteur du fichier contenant la variable */
  /*----------------------------------------------------------*/

  ECS_MALLOC(nom_fic_var,
             strlen(cas_ens->prefixe_rep) + strlen(cas_ens->prefixe_fic)
             + 1 + strlen(nom_var) + 1,
             char);

  sprintf(nom_fic_var, "%s%s.%s",
          cas_ens->prefixe_rep, cas_ens->prefixe_fic,
          nom_var);

  if (num_sor_prec == 0)
    mode_ouverture = ECS_FILE_MODE_APPEND;
  else
    mode_ouverture = ECS_FILE_MODE_WRITE;

  if (cas_ens->text == false) {

    fic_var = ecs_file_open(nom_fic_var,
                            mode_ouverture,
                            ECS_FILE_TYPE_BINARY);

    if (cas_ens->big_endian == true)
      ecs_file_set_big_endian(fic_var);

  }
  else
    fic_var = ecs_file_open(nom_fic_var,
                            mode_ouverture,
                            ECS_FILE_TYPE_TEXT);

  /* Info sur la création du fichier contenant la variable */
  /*-------------------------------------------------------*/

  printf("  %s %s\n", _("Creating file:"), nom_fic_var);

  ECS_FREE(nom_fic_var);

  /* Écriture dans ce fichier de la 1ère ligne (ligne de commentaires) */
  /*-------------------------------------------------------------------*/

  if (num_sor_prec != 0) {

    strcpy(buf, "Variable: ");
    strncpy(buf + strlen(buf), nom_var, 80 - strlen(buf));
    strncpy(buf + strlen(buf), ", Case: ", 80 - strlen(buf));
    strncpy(buf + strlen(buf), cas_ens->prefixe_fic, 80 - strlen(buf));

    ecs_post_ens__ecr_chaine(fic_var, buf);

  }

  ECS_FREE(nom_var);

  /* On réécrit le fichier Case */
  /*----------------------------*/

  if (cas_ens->modifie == true)
    ecs_loc_post_ens__ecr_case(cas_ens);

  return fic_var;
}

/*----------------------------------------------------------------------------*/


