#ifndef _ECS_POST_ENS_H_
#define _ECS_POST_ENS_H_

/*============================================================================
 *  Prototypes des fonctions de base
 *   réalisant les sorties pour post-traitement Ensight
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

/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie standard C ou BFT
 *----------------------------------------------------------------------------*/

#include <ecs_file.h>


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/

#include "ecs_def.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  des paquetages visibles
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 * Fichiers `include' publics  du  paquetage courant
 *----------------------------------------------------------------------------*/


/*============================================================================
 *                          Definitions de types
 *============================================================================*/

/* Structure définissant un cas EnSight */

typedef struct _ecs_post_ens_t ecs_post_ens_t;


/*============================================================================
 *                         Définitions d'énumérations
 *============================================================================*/


/*============================================================================
 *                       Définition de macros
 *============================================================================*/


/*============================================================================
 *                       Prototypes de fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction initialisant une structure `ecs_post_ens_t`.
 *----------------------------------------------------------------------------*/

ecs_post_ens_t  *
ecs_post_ens__cree_cas(const char  *nom_cas);

/*----------------------------------------------------------------------------
 *  Fonction détruisant une structure `ecs_post_ens_t`.
 *----------------------------------------------------------------------------*/

ecs_post_ens_t  *
ecs_post_ens__detruit_cas(ecs_post_ens_t  *cas_ens);

/*----------------------------------------------------------------------------
 *  Écriture d'une chaîne de caractères dans un fichier EnSight Gold.
 *----------------------------------------------------------------------------*/

void
ecs_post_ens__ecr_chaine(const ecs_file_t  *fic,
                         const char        *chaine);

/*----------------------------------------------------------------------------
 *  Écriture d'un entier dans un fichier EnSight Gold.
 *----------------------------------------------------------------------------*/

void
ecs_post_ens__ecr_int(const ecs_file_t  *fic,
                      int                val);

/*----------------------------------------------------------------------------
 *  Fonction écrivant le fichier contenant la géométrie
 *----------------------------------------------------------------------------*/

ecs_file_t  *
ecs_post_ens__ecrit_fic_geo(ecs_post_ens_t  *cas_ens);

/*----------------------------------------------------------------------------
 *  Fonction construisant le descripteur (`ecs_file_t') du fichier
 *  contenant les valeurs de la variable à sortir pour post-traitement Ensight
 *
 *  La fonction détermine aussi la ligne spécifiant la variable
 *   devant figurer dans le fichier Case
 *----------------------------------------------------------------------------*/

ecs_file_t  *
ecs_post_ens__ecrit_fic_var(ecs_post_ens_t  *cas_ens,
                            const char      *nom_var);

/*----------------------------------------------------------------------------*/

#endif /* _ECS_POST_ENS_H_ */
