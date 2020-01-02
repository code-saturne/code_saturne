#ifndef _ECS_PRE_GMSH_H_
#define _ECS_PRE_GMSH_H_

/*============================================================================
 *  Définition de la fonction
 *   de lecture d'un fichier de maillage au format Gmsh
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

/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie standard C ou BFT
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
 *  Lecture d'un fichier de maillage au format Gmsh
 *   et affectation des donnees dans la structure de maillage
 *----------------------------------------------------------------------------*/

ecs_maillage_t  *
ecs_pre_gmsh__lit_maillage(const char  *nom_fic_maillage);

/*----------------------------------------------------------------------------*/

#endif /* _ECS_PRE_GMSH_H_ */
