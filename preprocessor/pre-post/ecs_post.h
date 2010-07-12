#ifndef _ECS_POST_H_
#define _ECS_POST_H_

/*============================================================================
 *  Définition de types énumérés pour post-traitement
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
  ECS_POST_TYPE_INFO,     /* maillage d'information */
  ECS_POST_TYPE_ERREUR    /* maillage d'une zone avec erreur */

} ecs_post_type_t ;


/*============================================================================
 *                             Structures de données
 *============================================================================*/

/* Structure liée aux options amont de post-traitement */

typedef struct {

  bool     no_poly;          /* Suppression des polygones et des polyèdres */
  bool     text;             /* Forcer la sortie en mode texte si disponible */
  bool     big_endian;       /* Forcer la sortie binaire en mode big-endian
                                si possible */
  bool     color_to_group;   /* Forcer la sortie binaire en mode big-endian
                                si possible */
  bool     ecr_type[3];      /* Indicateur de sortie des maillages par type */

} ecs_post_opt_t;


/* Structure liée aux cas de post-traitement */

typedef struct {

  char            *nom_cas;           /* Nom du cas par défaut */

  bool             post_ens;          /* Indicateur post-traitement Ensight */
  ecs_post_opt_t   opt_ens;           /* Options pour EnSight */
  ecs_post_ens_t  *cas_ens;           /* Cas EnSight associé */

#if defined(HAVE_CGNS)

  bool             post_cgns;         /* Indicateur post-traitement CGNS */
  ecs_post_opt_t   opt_cgns;          /* Options pour CGNS */
  ecs_post_cgns_t *cas_cgns;          /* Cas CGNS associé */

#endif

#if defined(HAVE_MED)

  bool             post_med;          /* Indicateur post-traitement MED */
  ecs_post_opt_t   opt_med;           /* Options pour MED */
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
