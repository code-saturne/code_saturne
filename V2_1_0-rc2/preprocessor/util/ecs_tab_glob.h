#ifndef _ECS_TAB_GLOB_H_
#define _ECS_TAB_GLOB_H_

/*============================================================================
 *  Définition publique de la structure `tab_t' décrivant un tableau
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

#include <stdio.h>


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/

#include "ecs_def.h"


/*============================================================================
 *                         Définition d'enumération
 *============================================================================*/


/*============================================================================
 *                         Définitions de type
 *============================================================================*/

typedef struct {

  size_t        nbr;
  ecs_int_t    *val;

} ecs_tab_int_t;

typedef struct {

  size_t       nbr;
  char       **val;

} ecs_tab_char_t;


typedef struct {

  size_t       nbr;
  bool        *val;

} ecs_tab_bool_t;

/*----------------------------------------------------------------------------*/

#endif /* _ECS_TAB_GLOB_H_ */
