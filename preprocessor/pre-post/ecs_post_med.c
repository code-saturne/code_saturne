/*============================================================================
 *  Définitions des fonctions de base
 *   réalisant les sorties au format MED
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
ecs_post_med__cree_cas(const char  *nom_cas)
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

  lng_nom_fic = strlen(nom_cas) + strlen(".med") + 1;

  ECS_MALLOC(cas_med->nom_fic, lng_nom_fic, char);
  strcpy(cas_med->nom_fic, nom_cas);
  strcat(cas_med->nom_fic, ".med");

  for (ind = 0; ind < (ecs_int_t)strlen(nom_cas); ind++) {
    if (cas_med->nom_fic[ind] == ' ')
      cas_med->nom_fic[ind] = '_';
  }

  cas_med->nbr_maillages = 0;
  cas_med->tab_maillages = NULL;

  cas_med->fid = 0;

  for (ind = 0; ind < 3; ind++)  /* Useful for read only */
    cas_med->version[ind] = 0;

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

    if (cas_med->fid != 0) {
      if (MEDfileClose(cas_med->fid) != 0)
        ecs_error(__FILE__, __LINE__, 0,
                  _("MED: error closing file \"%s\"."),
                  cas_med->nom_fic);
    }

    ECS_FREE(cas_med->nom_cas);
    ECS_FREE(cas_med->nom_fic);

    for (ind = 0; ind < cas_med->nbr_maillages; ind++) {

      maillage_med = cas_med->tab_maillages[ind];

      ECS_FREE(maillage_med->nom_maillage);
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
  int         ind;
  ecs_int_t   lng_nom_maillage;

  ecs_med_maillage_t  *maillage_med;

  char  desc_maillage_med[MED_COMMENT_SIZE + 1] = "";

  char  dtunit[MED_LNAME_SIZE + 1] = "s";
  char  axisname[MED_SNAME_SIZE*3 + 1];
  char  axisunit[MED_SNAME_SIZE*3 + 1];

  med_err      ret_med = 0;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (cas_med == NULL)
    return;

  /* Create MED file if not done yet */

  if (cas_med->fid == 0) {

    cas_med->fid = MEDfileOpen(cas_med->nom_fic, MED_ACC_CREAT);

    if (cas_med->fid < 0)
      ecs_error(__FILE__, __LINE__, 0,
                _("MED: error opening file \"%s\"."),
                cas_med->nom_fic);

    printf("  %s %s\n", _("Creating file:"), cas_med->nom_fic);
  }

  /* Initialize strings */

  for (ind = 0; ind < MED_COMMENT_SIZE; ind++)
    desc_maillage_med[ind] = ' ';
  desc_maillage_med[MED_COMMENT_SIZE] = '\0';

  for (ind = 0; ind < MED_SNAME_SIZE*3; ind++) {
    axisname[ind] = ' ';
    axisunit[ind] = ' ';
  }
  axisname[0] = 'x';
  axisname[MED_SNAME_SIZE] = 'y';
  axisname[MED_SNAME_SIZE*2] = 'z';
  axisname[MED_SNAME_SIZE*3] = '\0';
  for (ind = 0; ind < 3; ind++)
    axisunit[ind * MED_SNAME_SIZE] = 'm';
  axisunit[MED_SNAME_SIZE*3] = '\0';

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

  strncpy(maillage_med->nom_maillage_med, nom_maillage, MED_NAME_SIZE);

  for (ind = 0; ind < lng_nom_maillage; ind++) {
    if (maillage_med->nom_maillage_med[ind] == ' ')
      maillage_med->nom_maillage_med[ind] = '_';
    else
      maillage_med->nom_maillage_med[ind]
        = tolower(maillage_med->nom_maillage_med[ind]);
  }
  for (ind = lng_nom_maillage; ind < MED_NAME_SIZE; ind++)
    maillage_med->nom_maillage_med[ind] = '\0';
  maillage_med->nom_maillage_med[MED_NAME_SIZE] = '\0';

  cas_med->nbr_maillages += 1;
  ECS_REALLOC(cas_med->tab_maillages,
              cas_med->nbr_maillages,
              ecs_med_maillage_t *);
  cas_med->tab_maillages[cas_med->nbr_maillages - 1] = maillage_med;

  /* Initialisation du maillage */

  ret_med = MEDmeshCr(cas_med->fid,
                      maillage_med->nom_maillage_med,
                      (med_int)3,
                      (med_int)(dim_m),
                      MED_UNSTRUCTURED_MESH,
                      desc_maillage_med,
                      dtunit,
                      MED_SORT_DTIT,
                      MED_CARTESIAN,
                      axisname,
                      axisunit);

  if (ret_med != 0)
    ecs_error(__FILE__, __LINE__, 0,
              _("MED: error writing file \"%s\".\n"
                "Name of mesh to create: \"%s\"\n"),
              cas_med->nom_fic, maillage_med->nom_maillage_med, (int)3);

  if (ret_med != 0)
    ecs_error(__FILE__, __LINE__, 0,
              _("MED: error writing file \"%s\".\n"
                "Mesh name: \"%s\"\n"
                "mesh dimension: \"%d\", space dimension: \"%d\""),
              cas_med->nom_fic, maillage_med->nom_maillage_med,
              (int)dim_m, 3);
}

/*----------------------------------------------------------------------------*/

#endif /* HAVE_MED */
