#ifndef _ECS_PRE_H_
#define _ECS_PRE_H_

/*============================================================================
 *  Definition des fonctions de base
 *   de remplissage de la structure de maillage a partir des donnees lues
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
 *                                 Visibilite
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie standard C
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdlib.h>
#include <string.h>


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/

#include "ecs_def.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles des paquetages visibles
 *----------------------------------------------------------------------------*/

#include "ecs_maillage.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage courant
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichier  `include' du  paquetage courant associe au fichier courant
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' prives   du  paquetage courant
 *----------------------------------------------------------------------------*/


/*============================================================================
 *                       Définition de type énuméré
 *============================================================================*/

typedef enum {

  ECS_PRE_FORMAT_NUL = -1,
  ECS_PRE_FORMAT_CGNS,
  ECS_PRE_FORMAT_CCM,
  ECS_PRE_FORMAT_COMET,
  ECS_PRE_FORMAT_ENS,
  ECS_PRE_FORMAT_GAMBIT,
  ECS_PRE_FORMAT_GMSH,
  ECS_PRE_FORMAT_IDEAS,
  ECS_PRE_FORMAT_MED,
  ECS_PRE_FORMAT_NOPO

} ecs_pre_format_t ;

/*============================================================================
 *                       Prototypes de fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction qui affiche la liste les formats supportés
 *----------------------------------------------------------------------------*/

void
ecs_pre__aff_formats(void);

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie le type de format de fichier associé à un fichier
 *   et à une clé optionnelle donnés
 *----------------------------------------------------------------------------*/

ecs_pre_format_t
ecs_pre__type_format(const char  *nom_fic,
                     const char  *mot_cle);

/*----------------------------------------------------------------------------
 *  Fonction qui lit les maillages sur fichiers
 *
 *  La fonction renvoie le maillage concaténé
 *----------------------------------------------------------------------------*/

ecs_maillage_t *
ecs_pre__lit_maillage(const char        *nom_fic,
                      ecs_pre_format_t   format,
                      int                num_maillage,
                      bool               cree_grp_cel_section,
                      bool               cree_grp_cel_zone,
                      bool               cree_grp_fac_section,
                      bool               cree_grp_fac_zone);

/*----------------------------------------------------------------------------*/

#endif /* _ECS_PRE_H_ */
