/*============================================================================
 *  Définitions des fonctions de base
 *   réalisant les sorties au format MED
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

#include "cs_config.h"

#if defined(HAVE_MED)


/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie standard C
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <ctype.h>
#include <string.h>


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/

#include "ecs_def.h"
#include "ecs_fic.h"
#include "ecs_mem.h"


/*----------------------------------------------------------------------------
 *  Fichier  `include' du  paquetage courant associe au fichier courant
 *----------------------------------------------------------------------------*/

#include "ecs_post_med.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' privés   du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_med_priv.h"



/*============================================================================
 *                              Fonctions privées
 *============================================================================*/


/*============================================================================
 *                             Fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction initialisant une structure `ecs_med_t` utilisée en écriture.
 *----------------------------------------------------------------------------*/

ecs_med_t *
ecs_post_med__cree_cas(const char  *nom_cas,
                       bool         no_poly)
{
  ecs_int_t    ind;
  ecs_int_t    lng_nom_fic;

  ecs_med_t  * cas_med;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Creation de la structure `ecs_med_t` et du fichier MED */
  /*--------------------------------------------------------*/

  ECS_MALLOC(cas_med, 1, ecs_med_t);

  ECS_MALLOC(cas_med->nom_cas, strlen(nom_cas) + 1, char);
  strcpy(cas_med->nom_cas, nom_cas);

  lng_nom_fic
    = strlen(nom_cas) + strlen(".med") + 1;

  ECS_MALLOC(cas_med->nom_fic, lng_nom_fic, char);
  strcpy(cas_med->nom_fic, nom_cas);
  strcat(cas_med->nom_fic, ".med");

  for (ind = 0; ind < (ecs_int_t)strlen(nom_cas); ind++) {
    if (cas_med->nom_fic[ind] == ' ')
      cas_med->nom_fic[ind] = '_';
  }

  cas_med->fid = MEDouvrir(cas_med->nom_fic, MED_CREATION);

  if (cas_med->fid < 0)
    ecs_error(__FILE__, __LINE__, 0,
              _("MED: error opening file \"%s\"."),
              cas_med->nom_fic);

  cas_med->nbr_var = 0;
  cas_med->nom_var = NULL;

  cas_med->nbr_maillages = 0;
  cas_med->tab_maillages = NULL;

  cas_med->no_poly = no_poly;

  /* Info sur la création du fichier MED */
  /*-------------------------------------*/

  printf("  %s %s\n", _("Creating file:"), cas_med->nom_fic);

  return cas_med;
}

/*----------------------------------------------------------------------------
 *  Fonction détruisant une structure `ecs_med_t` utilisée en écriture.
 *----------------------------------------------------------------------------*/

ecs_med_t *
ecs_post_med__detruit_cas(ecs_med_t  *cas_med)
{
  ecs_int_t ind;
  ecs_med_maillage_t  * maillage_med;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (cas_med != NULL) {

    /* Destruction du cas */

    if (MEDfermer(cas_med->fid) != 0)
      ecs_error(__FILE__, __LINE__, 0,
                _("MED: error closing file \"%s\"."),
                cas_med->nom_fic);

    ECS_FREE(cas_med->nom_cas);
    ECS_FREE(cas_med->nom_fic);

    for (ind = 0; ind < cas_med->nbr_var; ind++)
      ECS_FREE(cas_med->nom_var[ind]);

    ECS_FREE(cas_med->nom_var);
    cas_med->nbr_var = 0;

    for (ind = 0; ind < cas_med->nbr_maillages; ind++) {

      maillage_med = cas_med->tab_maillages[ind];

      ECS_FREE(maillage_med->nom_maillage);
      ECS_FREE(maillage_med->nbr_ele_typ);
      ECS_FREE(maillage_med->med_typ);

      ECS_FREE(maillage_med);

    }

    ECS_FREE(cas_med->tab_maillages);
    cas_med->nbr_maillages = 0;

    ECS_FREE(cas_med);

  }

  return cas_med;
}

/*----------------------------------------------------------------------------
 *  Fonction définissant un maillage pour une structure `ecs_med_t`.
 *----------------------------------------------------------------------------*/

void
ecs_post_med__ajoute_maillage(const char       *nom_maillage,
                              const ecs_int_t   dim_m,
                              ecs_med_t        *cas_med)
{
  ecs_int_t             ind;
  ecs_int_t             lng_nom_maillage;

  ecs_med_maillage_t  * maillage_med;

  char  desc_maillage_med[MED_TAILLE_DESC + 1] = "";

  med_err      ret_med = 0;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (cas_med == NULL)
    return;

  /* Vérification que le maillage n'a pas déjà été défini */

  for (ind = 0; ind < cas_med->nbr_maillages; ind++) {

    maillage_med = cas_med->tab_maillages[ind];

    if (strcmp(nom_maillage, maillage_med->nom_maillage) == 0)

      ecs_error(__FILE__, __LINE__, 0,
                _("A mesh named: %s\n"
                  "is already defined in MED case: %s\n"),
                nom_maillage, cas_med->nom_cas);

  }

  ECS_MALLOC(maillage_med, 1, ecs_med_maillage_t);

  lng_nom_maillage = strlen(nom_maillage);
  ECS_MALLOC(maillage_med->nom_maillage, lng_nom_maillage + 1, char);
  strcpy(maillage_med->nom_maillage, nom_maillage);

  strncpy(maillage_med->nom_maillage_med, nom_maillage, MED_TAILLE_NOM);

  for (ind = 0; ind < lng_nom_maillage; ind++) {
    if (maillage_med->nom_maillage_med[ind] == ' ')
      maillage_med->nom_maillage_med[ind] = '_';
    else
      maillage_med->nom_maillage_med[ind]
        = tolower(maillage_med->nom_maillage_med[ind]);
  }
  for (ind = lng_nom_maillage; ind < MED_TAILLE_NOM; ind++)
    maillage_med->nom_maillage_med[ind] = '\0';
  maillage_med->nom_maillage_med[MED_TAILLE_NOM] = '\0';

  /* BUG: En théorie, on devrait utiliser
     maillage_med->dim_entite = dim_entite;
     mais L'API MED est incohérente */

  maillage_med->dim_entite = 3;
  maillage_med->nbr_typ_ele = 0;
  maillage_med->nbr_ele_typ = NULL;
  maillage_med->med_typ = NULL;

  cas_med->nbr_maillages += 1;
  ECS_REALLOC(cas_med->tab_maillages,
              cas_med->nbr_maillages,
              ecs_med_maillage_t *);
  cas_med->tab_maillages[cas_med->nbr_maillages - 1] = maillage_med;

  /* Initialisation du maillage */

  ret_med = MEDmaaCr(cas_med->fid,
                     maillage_med->nom_maillage_med,
                     (med_int)(maillage_med->dim_entite),
                     MED_NON_STRUCTURE,
                     desc_maillage_med);

  if (ret_med != 0)
    ecs_error(__FILE__, __LINE__, 0,
              _("MED: error writing file \"%s\".\n"
                "Name of mesh to create: \"%s\"\n"),
              cas_med->nom_fic, maillage_med->nom_maillage_med, (int)3);

  ret_med = MEDdimEspaceCr(cas_med->fid,
                           maillage_med->nom_maillage_med,
                           (med_int)3);

  if (ret_med != 0)
    ecs_error(__FILE__, __LINE__, 0,
              _("MED: error writing file \"%s\".\n"
                "Mesh name: \"%s\"\n"
                "mesh dimension: \"%d\", space dimension: \"%d\""),
              cas_med->nom_fic, maillage_med->nom_maillage_med,
              (int)dim_m, 3);
}

/*----------------------------------------------------------------------------
 *  Fonction vérifiant si un champ figure dans la liste des champs du cas
 *----------------------------------------------------------------------------*/

bool
ecs_post_med__test_champ_liste(const char  *nom_champ,
                               ecs_med_t   *cas_med)
{
  ecs_int_t  ind;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for (ind = 0;
       ind < cas_med->nbr_var
         && strncmp(nom_champ, cas_med->nom_var[ind], MED_TAILLE_NOM) != 0;
       ind++);

  if (ind < cas_med->nbr_var)
    return true;
  else
    return false;
}

/*----------------------------------------------------------------------------
 *  Fonction ajoutant un champ à la liste des champs du cas
 *----------------------------------------------------------------------------*/

void
ecs_post_med__ajoute_champ_liste(const char  *nom_champ,
                                 ecs_med_t   *cas_med)
{
  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (ecs_post_med__test_champ_liste(nom_champ, cas_med) == true)
    ecs_error(__FILE__, __LINE__, 0,
              _("Field \"%s\" has already been defined in file \"%s\"\n"),
              nom_champ, cas_med->nom_fic);

  cas_med->nbr_var += 1;
  ECS_REALLOC(cas_med->nom_var, cas_med->nbr_var, char *);
  ECS_MALLOC(cas_med->nom_var[cas_med->nbr_var - 1], MED_TAILLE_NOM + 1, char);
  strncpy(cas_med->nom_var[cas_med->nbr_var - 1], nom_champ, MED_TAILLE_NOM);
  cas_med->nom_var[cas_med->nbr_var - 1][MED_TAILLE_NOM] = '\0';
}

/*----------------------------------------------------------------------------*/

#endif /* HAVE_MED */
