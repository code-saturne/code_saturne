#ifndef _ECS_FAMILLE_H_
#define _ECS_FAMILLE_H_

/*============================================================================
 *  Prototypes des fonctions de base
 *   associées à la structure `ecs_famille_t' décrivant une famille
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

/*---------------------------------------------------------------------------
 *  Fichiers `include' librairie standard C ou BFT
 *----------------------------------------------------------------------------*/

#include <stdio.h>


/*---------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/

#include "ecs_def.h"
#include "ecs_tab.h"


/*---------------------------------------------------------------------------
 *  Fichiers `include' publics  des paquetages visibles
 *----------------------------------------------------------------------------*/

#include "ecs_descr.h"


/*---------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage courant
 *----------------------------------------------------------------------------*/


/*============================================================================
 *                         Déclaration de la structure
 *============================================================================*/

typedef struct _ecs_famille_t ecs_famille_t;


/*============================================================================
 *                         Définition d'enumération
 *============================================================================*/

/*============================================================================
 *                       Prototypes de fonctions publiques
 *============================================================================*/

/*---------------------------------------------------------------------------
 *    Fonction de création d'une structure de famille `ecs_famille_t'
 *----------------------------------------------------------------------------*/

ecs_famille_t *
ecs_famille__cree(int           num,
                  ecs_descr_t  *descr_tete);

/*---------------------------------------------------------------------------
 *  Fonction libérant la structure `ecs_famille_t' donnée en argument.
 *  Elle renvoie un pointeur NULL
 *----------------------------------------------------------------------------*/

ecs_famille_t *
ecs_famille__detruit(ecs_famille_t  *this_fam);

/*---------------------------------------------------------------------------
 *  Fonction imprimant le contenu d'une structure `ecs_famille_t' donnée
 *   sur le flux décrit par la structure `ecs_file_t'
 *----------------------------------------------------------------------------*/

void
ecs_famille__imprime(const ecs_famille_t  *this_fam,
                     int                   imp_col,
                     FILE                 *fic_imp);

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie la taille en octets d'une structure `ecs_famille_t'
 *----------------------------------------------------------------------------*/

float
ecs_famille__ret_taille(const ecs_famille_t  *this_fam);

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie le nombre et la liste des des pointeurs sur les noms
 *   des descripteurs de la famille donnée en argument
 *----------------------------------------------------------------------------*/

ecs_tab_char_t
ecs_famille__ret_nom(const ecs_famille_t  *this_fam);

/*----------------------------------------------------------------------------
 *  Fonction qui alloue une structure `ecs_famille_t' et qui remplit
 *   son contenu en copiant le contenu de la structure donnée en argument
 *----------------------------------------------------------------------------*/

ecs_famille_t *
ecs_famille__copie(ecs_famille_t  *this_famille);

/*----------------------------------------------------------------------------
 *  Fonction qui affiche la définition de la famille
 *----------------------------------------------------------------------------*/

void
ecs_famille__affiche(const ecs_famille_t  *this_fam);

/*----------------------------------------------------------------------------*/

#endif /* _ECS_FAMILLE_H_ */
