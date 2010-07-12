#ifndef _ECS_DESCR_PRIV_H_
#define _ECS_DESCR_PRIV_H_

/*============================================================================
 *  Définition privee de la structure `ecs_descr_t' décrivant
 *   un descripteur de champ de type "attribut"
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


/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie standard C
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/

#include "ecs_def.h"
#include "ecs_tab_glob.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  des paquetages visibles
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_descr.h"


/*============================================================================
 *                       Définition de macros
 *============================================================================*/


/*============================================================================
 *                       Définition de la structure
 *============================================================================*/

/* Descripteur de couleur ou groupe. */

struct _ecs_descr_t {

  int                  num;           /* Numéro */
  ecs_descr_typ_t      typ;           /* Type */
  int                  ide;           /* Numéro d'attribut */
  char               * nom;           /* Nom */
  struct
  _ecs_descr_t       * l_descr_sui;   /* Pointeur sur le descripteur suivant */

};

/*----------------------------------------------------------------------------*/

#endif /* _ECS_DESCR_PRIV_H_ */
