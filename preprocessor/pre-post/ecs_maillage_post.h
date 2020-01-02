#ifndef _ECS_MAILLAGE_POST_H_
#define _ECS_MAILLAGE_POST_H_

/*============================================================================
 *  Prototypes des fonctions
 *   associées à la structure `ecs_maillage_t' décrivant un maillage
 *   et réalisant les sorties pour post-traitement
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
 *                                 Visibilite
 *============================================================================*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie standard C
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/

#include "ecs_tab.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage global "Post-Traitement"
 *----------------------------------------------------------------------------*/

#include "ecs_post.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  des paquetages visibles
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_maillage.h"


/*============================================================================
 *                       Prototypes de fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *    Fonction qui écrit le maillage sur fichier pour post-traitement
 *  Seuls les éléments de l'entité de rang `ent_post' sont écrits.
 *
 *    Les éléments de l'entité à écrire doivent être sous forme de
 *  connectivité nodale.
 *----------------------------------------------------------------------------*/

void
ecs_maillage_post__ecr(const char       *nom_maillage,
                       ecs_maillage_t   *maillage,
                       ecs_post_type_t   type_post,
                       ecs_post_t       *cas_post);

/*----------------------------------------------------------------------------
 *  Écriture du maillage correspondant à une liste de faces sur fichier
 *  pour post-traitement.
 *
 *  Cette fonction crée une coupe correspondant à la liste de faces donnée
 *  (ce qui provoque automatiquement son post-traitement), puis la détruit.
 *  Le nom utilisé pour cette sortie ne sera donc plus disponible pour
 *  d'autres coupes.
 *
 *  Le maillage principal doit être en connectivité descendante.
 *----------------------------------------------------------------------------*/

void
ecs_maillage_post__ecr_fac_liste(const char           *nom_liste,
                                 ecs_maillage_t       *maillage,
                                 const ecs_tab_int_t   liste_fac,
                                 ecs_post_type_t       type_post,
                                 ecs_post_t           *cas_post);

/*----------------------------------------------------------------------------*/

#endif /* _ECS_MAILLAGE_POST_H_ */
