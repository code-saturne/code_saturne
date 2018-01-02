#ifndef _ECS_POST_H_
#define _ECS_POST_H_

/*============================================================================
 *  Définition de types énumérés pour post-traitement
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


/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie standard C
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/

#include "ecs_def.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  des paquetages visibles
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_post_ens.h"
#include "ecs_med.h"
#include "ecs_post_cgns.h"


/*============================================================================
 *                         Définitions d'énumerations
 *============================================================================*/

/*  Définition d'énumération liée au type de post-traitement */

typedef enum {

  ECS_POST_TYPE_VOLUME,   /* maillage volumique (initial) */
  ECS_POST_TYPE_ERREUR    /* maillage d'une zone avec erreur */

} ecs_post_type_t ;


/*============================================================================
 *                             Structures de données
 *============================================================================*/

/* Structure liée aux cas de post-traitement */

typedef struct {

  char            *nom_cas;           /* Nom du cas par défaut */

  bool             opt_ens[2];        /* Indicateur de sortie par type */
  ecs_post_ens_t  *cas_ens;           /* Cas EnSight associé */

#if defined(HAVE_CGNS)

  bool             opt_cgns[2];       /* Indicateur de sortie par type */
  ecs_post_cgns_t *cas_cgns;          /* Cas CGNS associé */

#endif

#if defined(HAVE_MED)

  bool             opt_med[2];        /* Indicateur de sortie par type */
  ecs_med_t       *cas_med;           /* Cas MED associé */

#endif

} ecs_post_t;

/*============================================================================
 *                       Prototypes de fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction initialisant une structure `ecs_post_t`
 *----------------------------------------------------------------------------*/

ecs_post_t *
ecs_post__cree_cas(const char  *nom_cas);

/*----------------------------------------------------------------------------
 *  Fonction détruisant une structure `ecs_post_t`
 *----------------------------------------------------------------------------*/

ecs_post_t *
ecs_post__detruit_cas(ecs_post_t  *cas);

/*----------------------------------------------------------------------------*/

#endif /* _ECS_POST_H_ */
