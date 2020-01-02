#ifndef _ECS_POST_CGNS_H_
#define _ECS_POST_CGNS_H_

/*============================================================================
 *  Prototypes des fonctions de base
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


/*============================================================================
 *                         Définitions d'énumerations
 *============================================================================*/


/*============================================================================
 * Définitions de types
 *============================================================================*/

/* Structure définissant un cas CGNS */

typedef struct _ecs_post_cgns_t ecs_post_cgns_t;

/*============================================================================
 *                           Définition de macros
 *============================================================================*/


/*============================================================================
 *                       Prototypes de fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction initialisant une structure `ecs_cgns_t` utilisée en ecriture.
 *---------------------------------------------------------------------------*/

ecs_post_cgns_t  *
ecs_post_cgns__cree_cas(const char  *nom_cas);

/*----------------------------------------------------------------------------
 *  Fonction détruisant une structure `ecs_cgns_t` utilisée en écriture.
 *----------------------------------------------------------------------------*/

ecs_post_cgns_t  *
ecs_post_cgns__detruit_cas(ecs_post_cgns_t  *cas_cgns);

/*----------------------------------------------------------------------------
 *  Fonction fermant le fichier associé à un cas CGNS
 *  (pour forcer sa mise à jour).
 *---------------------------------------------------------------------------*/

void
ecs_post_cgns__ferme_cas(ecs_post_cgns_t  *cas_cgns);

/*----------------------------------------------------------------------------
 *  Fonction définissant un maillage pour une structure `ecs_cgns_t`.
 *----------------------------------------------------------------------------*/

void
ecs_post_cgns__ajoute_maillage(const char       *nom_maillage,
                               ecs_int_t         dim_entite,
                               ecs_post_cgns_t  *cas_cgns);

#endif /* HAVE_CGNS */

/*----------------------------------------------------------------------------*/

#endif /* _ECS_POST_CGNS_H_ */
