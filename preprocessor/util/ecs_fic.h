#ifndef _ECS_FIC_H_
#define _ECS_FIC_H_

/*============================================================================
 *  Prototypes des fonctions
 *   associées aux impressions dans un fichier
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

#include "ecs_def.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction d'impression du nom et de la valeur d'un pointeur
 *----------------------------------------------------------------------------*/

void
ecs_fic__imprime_ptr(FILE        *fic_imp,
                     int          profondeur_imp,
                     const char  *nom_ptr,
                     const void  *ptr);

/*----------------------------------------------------------------------------
 *  Fonction d'impression du nom et de la valeur d'une variable
 *----------------------------------------------------------------------------*/

void
ecs_fic__imprime_val(FILE        *fic_imp,
                     int          profondeur_imp_nom,
                     const char  *nom,
                     ecs_type_t   typ_e,
                     const void  *val);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* _ECS_FIC_H_ */
