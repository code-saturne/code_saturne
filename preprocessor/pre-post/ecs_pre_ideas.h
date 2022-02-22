#ifndef _ECS_PRE_IDEAS_H_
#define _ECS_PRE_IDEAS_H_

/*============================================================================
 *  Définition de la fonction
 *   de lecture d'un fichier de maillage IDEAS-MS au format "universel"
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


/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie standard C
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles des paquetages visibles
 *----------------------------------------------------------------------------*/

#include "ecs_maillage.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage courant
 *----------------------------------------------------------------------------*/


/*============================================================================
 *  Prototypes de fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Lecture d'un fichier I-DEAS Master Series au format universel
 *   et affectation des donnees dans la structure de maillage
 *
 *  Hypothèses de lecture :
 *   le dataset sur les groupes doit se trouver placer dans le fichier Ideas
 *   après les dataset sur les noeuds et sur les éléments ; le dataset
 *   sur noeuds doit lui-même se trouver après celui sur les systèmes de
 *   coordonnées.
*----------------------------------------------------------------------------*/

ecs_maillage_t *
ecs_pre_ideas__lit_maillage(const char  *nom_fic_maillage);

/*----------------------------------------------------------------------------*/

#endif /* _ECS_PRE_IDEAS_H_ */
