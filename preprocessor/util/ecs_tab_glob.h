#ifndef _ECS_TAB_GLOB_H_
#define _ECS_TAB_GLOB_H_

/*============================================================================
 *  Définition publique de la structure `tab_t' décrivant un tableau
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2011 EDF S.A.

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
