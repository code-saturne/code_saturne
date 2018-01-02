#ifndef _ECS_DESCR_H_
#define _ECS_DESCR_H_

/*============================================================================
 *  Prototypes des fonctions de base
 *   associées à la structure `ecs_descr_t' décrivant un descripteur de table
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
 *                                 Visibilité
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie standard C ou BFT
 *----------------------------------------------------------------------------*/

#include <stdio.h>


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/

#include "ecs_def.h"
#include "ecs_tab.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  des paquetages visibles
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage courant
 *----------------------------------------------------------------------------*/


/*============================================================================
 *                         Déclaration de la structure
 *============================================================================*/

typedef struct _ecs_descr_t ecs_descr_t;


/*============================================================================
 *                         Définition de macros
 *============================================================================*/

#define ECS_DESCR_NUM_NUL      -1

#define ECS_DESCR_ID_EXT_NUL   -1
#define ECS_DESCR_IDE_NUL      -1


/*============================================================================
 *                       Prototypes de fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction de création d'une structure de descripteur de table
 *----------------------------------------------------------------------------*/

ecs_descr_t *
ecs_descr__cree(int          ide,
                const char  *nom);

/*----------------------------------------------------------------------------
 *  Fonction libérant la structure `ecs_descr_t' donnée en argument.
 *  Elle renvoie un pointeur NULL
 *----------------------------------------------------------------------------*/

ecs_descr_t *
ecs_descr__detruit(ecs_descr_t  *this_descr);

/*----------------------------------------------------------------------------
 *  Fonction imprimant le contenu d'une structure `ecs_descr_t' donnée
 *   sur le flux décrit par la structure `ecs_file_t'
 *----------------------------------------------------------------------------*/

void
ecs_descr__imprime(const ecs_descr_t  *this_descr,
                   ecs_int_t           imp_col,
                   FILE               *fic_imp);

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie la taille en octets d'une structure `ecs_descr_t'
 *----------------------------------------------------------------------------*/

float
ecs_descr__ret_taille(const ecs_descr_t *this_descr);

/*----------------------------------------------------------------------------
 *  Fonction qui alloue une structure `ecs_descr_t' et qui remplit
 *   son contenu en copiant le contenu de la structure donnée en argument
 *----------------------------------------------------------------------------*/

ecs_descr_t *
ecs_descr__copie(ecs_descr_t  *this_descr);

/*----------------------------------------------------------------------------
 *  Fonction qui compare 2 descripteurs
 *
 *  La fonction renvoie :
 *  - `true'  si les deux descripteurs sont identiques
 *  - `false' sinon
 *----------------------------------------------------------------------------*/

bool
ecs_descr__compare(const ecs_descr_t  *descr_1,
                   const ecs_descr_t  *descr_2);

/*----------------------------------------------------------------------------
 *  Fonction qui affiche le nom d'un descripteur
 *----------------------------------------------------------------------------*/

void
ecs_descr__affiche(const ecs_descr_t  *descr,
                   int                 decal);

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie le nom du descripteur donné en argument
 *----------------------------------------------------------------------------*/

const char *
ecs_descr__ret_nom(const ecs_descr_t  *descr);

/*----------------------------------------------------------------------------*/

#endif /* _ECS_DESCR_H_ */
