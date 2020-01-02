#ifndef _ECS_PRE_CGNS_H_
#define _ECS_PRE_CGNS_H_

/*============================================================================
 *  Définition de la fonction
 *   de lecture d'un fichier de maillage au format CGNS
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
 *  Fichiers `include' visibles du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/

#include "ecs_def.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage global "CGNS"
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles des paquetages visibles
 *----------------------------------------------------------------------------*/

#include "ecs_maillage.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage courant
 *----------------------------------------------------------------------------*/


/*============================================================================
 *                       Prototypes de fonctions publiques
 *============================================================================*/


/*----------------------------------------------------------------------------
 *  Lecture d'un fichier au format CGNS
 *   et affectation des données dans la structure de maillage
 *----------------------------------------------------------------------------*/

ecs_maillage_t *
ecs_pre_cgns__lit_maillage(const char  *nom_fic_maillage,
                           int          num_maillage,
                           bool         cree_grp_cel_section,
                           bool         cree_grp_cel_zone,
                           bool         cree_grp_fac_section,
                           bool         cree_grp_fac_zone);

#endif /* HAVE_CGNS */

/*----------------------------------------------------------------------------*/

#endif /* _ECS_PRE_CGNS_H_ */
