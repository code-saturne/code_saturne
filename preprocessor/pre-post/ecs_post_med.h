#ifndef _ECS_POST_MED_H_
#define _ECS_POST_MED_H_

/*============================================================================
 *  Prototypes des fonctions de base
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

#include "ecs_med.h"


/*============================================================================
 *                         Définitions d'énumerations
 *============================================================================*/


/*============================================================================
 *                       Définition de macros
 *============================================================================*/


/*============================================================================
 *                       Prototypes de fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction initialisant une structure `ecs_med_t` utilisée en écriture.
 *----------------------------------------------------------------------------*/

ecs_med_t *
ecs_post_med__cree_cas(const char  *nom_cas,
                       bool         no_poly);

/*----------------------------------------------------------------------------
 *  Fonction détruisant une structure `ecs_med_t` utilisée en écriture.
 *----------------------------------------------------------------------------*/

ecs_med_t *
ecs_post_med__detruit_cas(ecs_med_t  *cas_med);

/*----------------------------------------------------------------------------
 *  Fonction définissant un maillage pour une structure `ecs_med_t`.
 *----------------------------------------------------------------------------*/

void
ecs_post_med__ajoute_maillage(const char       *nom_maillage,
                              const ecs_int_t   dim_m,
                              ecs_med_t        *cas_med);


/*----------------------------------------------------------------------------
 *  Fonction vérifiant si un champ figure dans la liste des champs du cas
 *----------------------------------------------------------------------------*/

bool
ecs_post_med__test_champ_liste(const char  *nom_champ,
                               ecs_med_t   *cas_med);


/*----------------------------------------------------------------------------
 *  Fonction ajoutant un champ à la liste des champs du cas
 *----------------------------------------------------------------------------*/

void
ecs_post_med__ajoute_champ_liste(const char  *nom_champ,
                                 ecs_med_t   *cas_med);

#endif /* HAVE_MED */

/*----------------------------------------------------------------------------*/

#endif /* _ECS_POST_MED_H_ */

