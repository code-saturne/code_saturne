#ifndef _ECS_FAMILLE_PRIV_H_
#define _ECS_FAMILLE_PRIV_H_

/*============================================================================
 *  Définition privée de la structure `ecs_famille_t' décrivant une famille
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


/*============================================================================
 *                                 Visibilité
 *============================================================================*/

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

#include "ecs_descr.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_famille.h"


/*============================================================================
 *                       Définition de macros
 *============================================================================*/


/*============================================================================
 *                       Définition de la structure
 *============================================================================*/

struct _ecs_famille_t {

  int              num;            /* Numéro de la famille */
  ecs_descr_t     *descr;          /* Tête de la liste chaînee de
                                      descripteurs définissant la famille */
  struct
  _ecs_famille_t  *l_famille_sui;  /* Pointeur sur la famille suivante de la
                                      liste chaînée des familles */

};

/*============================================================================
 *  Les valeurs des tableaux `famille' sont les numéros des familles
 *   numerotées à partir de `1' comme suit :
 *   - en commencant par les familles des cellules
 *   - puis          par les familles des faces
 *============================================================================*/

/*----------------------------------------------------------------------------*/

#endif /* _ECS_FAMILLE_PRIV_H_ */
