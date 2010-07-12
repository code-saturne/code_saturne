#ifndef _ECS_CMD_H_
#define _ECS_CMD_H_

/*============================================================================
 *  Prototypes des fonctions de base
 *   associées à la structure `ecs_cmd_t' décrivant
 *   les options de la ligne de commande
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


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  des paquetages visibles
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage courant
 *----------------------------------------------------------------------------*/


/*============================================================================
 *                         Déclaration de la structure
 *============================================================================*/

typedef struct _ecs_cmd_post_t    ecs_cmd_post_t;
typedef struct _ecs_cmd_t         ecs_cmd_t;

/*============================================================================
 *                       Prototypes de fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction qui lit la ligne de commande
 *----------------------------------------------------------------------------*/

ecs_cmd_t *
ecs_cmd__lit_arg(int    argc,
                 char  *argv[]);

/*----------------------------------------------------------------------------
 *  Fonction liberant une structure `ecs_cmd_t' donnee en argument.
 *  Elle renvoie un pointeur NULL
 *----------------------------------------------------------------------------*/

ecs_cmd_t *ecs_cmd__detruit(ecs_cmd_t  *cmd);

/*----------------------------------------------------------------------------*/

#endif /* _ECS_CMD_H_ */
