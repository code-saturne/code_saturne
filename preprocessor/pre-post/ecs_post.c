/*============================================================================
 *  Définitions des fonctions de base
 *   réalisant les sorties pour post-traitement
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

#include "cs_config.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie standard C
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <string.h>


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/

#include "ecs_def.h"
#include "ecs_mem.h"


/*----------------------------------------------------------------------------
 *  Fichier  `include' du  paquetage courant associe au fichier courant
 *----------------------------------------------------------------------------*/

#include "ecs_post.h"
#include "ecs_post_ens.h"
#include "ecs_post_cgns.h"
#include "ecs_post_med.h"


/*============================================================================
 *                              Fonctions privées
 *============================================================================*/

/*============================================================================
 *                             Fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction initialisant une structure `ecs_post_t`
 *----------------------------------------------------------------------------*/

ecs_post_t *
ecs_post__cree_cas(const char  *nom_cas)
{
  ecs_post_t *cas;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(nom_cas != NULL);

  ECS_MALLOC(cas, 1, ecs_post_t);

  ECS_MALLOC(cas->nom_cas, strlen(nom_cas) + 1, char);
  strcpy(cas->nom_cas, nom_cas);

  cas->cas_ens = NULL;
  cas->opt_ens[ECS_POST_TYPE_VOLUME] = false;
  cas->opt_ens[ECS_POST_TYPE_ERREUR] = false;

#if defined(HAVE_CGNS)

  cas->cas_cgns = NULL;
  cas->opt_cgns[ECS_POST_TYPE_VOLUME] = false;
  cas->opt_cgns[ECS_POST_TYPE_ERREUR] = false;

#endif /* HAVE_CGNS */

#if defined(HAVE_MED)

  cas->cas_med = NULL;
  cas->opt_med[ECS_POST_TYPE_VOLUME] = false;
  cas->opt_med[ECS_POST_TYPE_ERREUR] = false;

#endif /* HAVE_MED */

  return cas;
}

/*----------------------------------------------------------------------------
 *  Fonction détruisant une structure `ecs_post_t`
 *----------------------------------------------------------------------------*/

ecs_post_t *
ecs_post__detruit_cas(ecs_post_t  *cas)
{
  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (cas == NULL)
    return NULL;

  if (cas->cas_ens != NULL)
    ecs_post_ens__detruit_cas(cas->cas_ens);

#if defined(HAVE_CGNS)

  if (cas->cas_cgns != NULL)
    ecs_post_cgns__detruit_cas(cas->cas_cgns);

#endif

#if defined(HAVE_MED)

  if (cas->cas_med != NULL)
    ecs_post_med__detruit_cas(cas->cas_med);

#endif

  ECS_FREE(cas->nom_cas);

  ECS_FREE(cas);

  return NULL;
}

/*----------------------------------------------------------------------------*/


