#ifndef _ECS_TABLE_COMM_H_
#define _ECS_TABLE_COMM_H_

/*============================================================================
 *  Prototypes des fonctions
 *   associées à la structure `ecs_table_t' décrivant une table
 *   et réalisant des communications
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


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage global "Post-Traitement"
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  des paquetages visibles
 *----------------------------------------------------------------------------*/

#include "ecs_comm.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_table.h"


/*============================================================================
 *                       Prototypes de fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction qui ecrit le tableau des positions d'une table
 *   dans le fichier d'interface pour le Noyau
 *----------------------------------------------------------------------------*/

void
ecs_table_comm__ecr_pos(ecs_table_t  *this_table,
                        const char   *comm_nom_rubrique,
                        size_t        location_id,
                        size_t        index_id,
                        ecs_comm_t   *comm);

/*----------------------------------------------------------------------------
 *  Fonction qui écrit le contenu d'une table
 *   dans le fichier d'interface pour le Noyau
 *----------------------------------------------------------------------------*/

void ecs_table_comm__ecr(ecs_table_t  *this_table,
                         const char   *comm_nom_rubrique,
                         size_t        location_id,
                         size_t        index_id,
                         size_t        n_location_values,
                         ecs_comm_t   *comm);

/*----------------------------------------------------------------------------*/

#endif /* _ECS_TABLE_COMM_H_ */
