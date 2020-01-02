/*============================================================================
 *  Définitions des fonctions de base
 *   réalisant les sorties au format CGNS
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

#include "cs_config.h"

#if defined(HAVE_CGNS)


/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie standard C
 *---------------------------------------------------------------------------*/

#include <assert.h>
#include <ctype.h>
#include <string.h>


/*----------------------------------------------------------------------------
 * Fichiers `include' publics  du  paquetage global "CGNS"
 *----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif

#include <cgnslib.h>

#ifdef __cplusplus
}
#endif


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage global "Utilitaire"
 *---------------------------------------------------------------------------*/

#include "ecs_def.h"
#include "ecs_fic.h"
#include "ecs_mem.h"


/*----------------------------------------------------------------------------
  *  Fichier  `include' du  paquetage courant associe au fichier courant
 *---------------------------------------------------------------------------*/

#include "ecs_post_cgns.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' privés   du  paquetage courant
 *---------------------------------------------------------------------------*/

#include "ecs_post_cgns_priv.h"


/*============================================================================
 *                  Définitions de paramètres et macros
 *============================================================================*/


/*============================================================================
 *                              Fonctions privées
 *============================================================================*/


/*============================================================================
 *                             Fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction initialisant une structure `ecs_cgns_t` utilisée en ecriture.
 *---------------------------------------------------------------------------*/

ecs_post_cgns_t  *
ecs_post_cgns__cree_cas(const char  *nom_cas)
{
  ecs_post_cgns_t  *cas_cgns;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Creation de la structure `ecs_cgns_t` */
  /*---------------------------------------*/

  ECS_MALLOC(cas_cgns, 1, ecs_post_cgns_t);

  ECS_MALLOC(cas_cgns->nom_cas, strlen(nom_cas) + 1, char);
  strcpy(cas_cgns->nom_cas, nom_cas);

  cas_cgns->nbr_bases = 0;
  cas_cgns->tab_bases = NULL;

  return cas_cgns;
}

/*----------------------------------------------------------------------------
 *  Fonction détruisant une structure `ecs_post_cgns_t` utilisée en écriture.
 *---------------------------------------------------------------------------*/

ecs_post_cgns_t  *
ecs_post_cgns__detruit_cas(ecs_post_cgns_t  *cas_cgns)
{
  ecs_int_t  ind;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (cas_cgns != NULL) {

    ECS_FREE(cas_cgns->nom_cas);

    /* Destruction des bases */

    for (ind = 0; ind < cas_cgns->nbr_bases; ind++) {

      if (cas_cgns->tab_bases[ind]->fic_ouvert == true) {

        if (cg_close(cas_cgns->tab_bases[ind]->num_fic) != CG_OK)
          ecs_error(__FILE__, __LINE__, 0,
                    _("CGNS: error closing file \"%s\"\n%s"),
                    cas_cgns->tab_bases[ind]->nom_fic, cg_get_error());

      }

      ECS_FREE(cas_cgns->tab_bases[ind]->nom_fic);
      ECS_FREE(cas_cgns->tab_bases[ind]->nom_maillage);

      ECS_FREE(cas_cgns->tab_bases[ind]);

    }

    ECS_FREE(cas_cgns->tab_bases);
    cas_cgns->nbr_bases = 0;

    ECS_FREE(cas_cgns);

  }

  return cas_cgns;
}

/*----------------------------------------------------------------------------
 *  Fonction fermant le fichier associé à un cas CGNS
 *  (pour forcer sa mise à jour).
 *---------------------------------------------------------------------------*/

void
ecs_post_cgns__ferme_cas(ecs_post_cgns_t  *cas_cgns)
{
  ecs_int_t  ind;
  ecs_post_cgns_base_t  * base_cgns;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for (ind = 0; ind < cas_cgns->nbr_bases; ind++) {

    base_cgns = cas_cgns->tab_bases[ind];

    if (base_cgns->fic_ouvert == true) {

      if (cg_close(base_cgns->num_fic) != CG_OK)
        ecs_error(__FILE__, __LINE__, 0,
                  _("CGNS: error closing file \"%s\"\n%s"),
                  base_cgns->nom_fic, cg_get_error());

      base_cgns->fic_ouvert = false;

    }
  }
}

/*----------------------------------------------------------------------------
 *  Fonction définissant un maillage pour une structure `ecs_post_cgns_t`.
 *----------------------------------------------------------------------------*/

void
ecs_post_cgns__ajoute_maillage(const char       *nom_maillage,
                               ecs_int_t         dim_entite,
                               ecs_post_cgns_t  *cas_cgns)
{
  ecs_int_t          ind;
  ecs_int_t          lng_nom_fic;
  ecs_int_t          lng_nom_maillage;

  ecs_post_cgns_base_t  * base_cgns;

  char               nom_maillage_cgns[ECS_CGNS_TAILLE_NOM + 1];

  int  num_base = 0;
  int  ret_cgns = 0;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (cas_cgns == NULL)
    return;

  /* Vérification que le maillage n'a pas déjà été défini */

  for (ind = 0; ind < cas_cgns->nbr_bases; ind++) {

    base_cgns = cas_cgns->tab_bases[ind];

    if (strcmp(nom_maillage, base_cgns->nom_maillage) == 0)

      ecs_error(__FILE__, __LINE__, 0,
                _("A mesh named: %s\n"
                  "is already defined in CGNS case: %s\n"),
                nom_maillage, cas_cgns->nom_cas);
  }

  /* Initialisation du maillage */
  /*----------------------------*/

  ECS_MALLOC(base_cgns, 1, ecs_post_cgns_base_t);

  lng_nom_maillage = strlen(nom_maillage);
  ECS_MALLOC(base_cgns->nom_maillage, lng_nom_maillage + 1, char);
  strcpy(base_cgns->nom_maillage, nom_maillage);

  strncpy(nom_maillage_cgns, nom_maillage, ECS_CGNS_TAILLE_NOM);
  nom_maillage_cgns[ECS_CGNS_TAILLE_NOM] = '\0';

  for (ind = lng_nom_maillage; ind < ECS_CGNS_TAILLE_NOM; ind++)
    nom_maillage_cgns[ind] = '\0';
  nom_maillage_cgns[ECS_CGNS_TAILLE_NOM] = '\0';

  /* Dimension espace et maillage */

  base_cgns->dim_entite = dim_entite;
  base_cgns->dim_espace = 3;

  /* Champs associés aux sorties */

  cas_cgns->nbr_bases += 1;
  ECS_REALLOC(cas_cgns->tab_bases,
              cas_cgns->nbr_bases,
              ecs_post_cgns_base_t *);
  cas_cgns->tab_bases[cas_cgns->nbr_bases - 1] = base_cgns;

  /* Création du fichier CGNS associé */
  /*----------------------------------*/

  lng_nom_fic =   strlen(cas_cgns->nom_cas) + 1
                + lng_nom_maillage + strlen(".cgns") + 1;

  ECS_MALLOC(base_cgns->nom_fic, lng_nom_fic, char);
  sprintf(base_cgns->nom_fic, "%s_%s.cgns", cas_cgns->nom_cas, nom_maillage);

  for (ind = 0; ind < lng_nom_fic; ind++) {
    if (base_cgns->nom_fic[ind] == ' ')
      base_cgns->nom_fic[ind] = '_';
    else
      base_cgns->nom_fic[ind] = tolower(base_cgns->nom_fic[ind]);
  }

  if (   cg_open(base_cgns->nom_fic, CG_MODE_WRITE, &(base_cgns->num_fic))
      != CG_OK)
    ecs_error(__FILE__, __LINE__, 0,
              _("CGNS: error opening file \"%s\":\n%s"),
              base_cgns->nom_fic, cg_get_error());

  /* Info sur la creation du fichier CGNS */
  /*-------------------------------------*/

  printf("  %s %s\n", _("Creating file:"), base_cgns->nom_fic);

  base_cgns->fic_ouvert = true;

  /* Ajout d'une base associée au cas CGNS */

  ret_cgns = cg_base_write(base_cgns->num_fic, nom_maillage_cgns,
                           (int)dim_entite, 3, &num_base);

  if (ret_cgns != CG_OK)
    ecs_error(__FILE__, __LINE__, 0,
              _("CGNS: error writing file \"%s\":\n"
                "Name of mesh to write: \"%s\"\n%s"),
              base_cgns->nom_fic, nom_maillage_cgns, cg_get_error());

  if (num_base != 1)
    ecs_error(__FILE__, __LINE__, 0,
              _("CGNS: error writing file \"%s\":\n"
                "The returned base number is \"%d\" and not \"1\"."),
              base_cgns->nom_fic, num_base);
}

/*----------------------------------------------------------------------------*/

#endif /* HAVE_CGNS */
