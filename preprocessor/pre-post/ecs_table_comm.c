/*============================================================================
 *  Definitions des fonctions
 *   associees a la structure `ecs_table_t' decrivant une table
 *   et realisant les sorties pour la communication
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
 *                                 Visibilite
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie standard C
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <string.h> /* strlen() */


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/

#include "ecs_def.h"
#include "ecs_mem.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage global "Post-Traitement"
 *----------------------------------------------------------------------------*/

#include "ecs_comm.h"


/*----------------------------------------------------------------------------
  *  Fichiers `include' visibles des paquetages visibles
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_table.h"


/*----------------------------------------------------------------------------
 *  Fichier  `include' du  paquetage courant associe au fichier courant
 *----------------------------------------------------------------------------*/

#include "ecs_table_comm.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' prives   du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_table_priv.h"


/*============================================================================
 *                              Fonctions privees
 *============================================================================*/

/*============================================================================
 *                             Fonctions publiques
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
                        ecs_comm_t   *comm)
{
  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(this_table != NULL);

  ecs_table__regle_en_pos(this_table);

  ecs_comm_write_section(comm_nom_rubrique,
                         this_table->nbr + 1,
                         location_id,
                         index_id,
                         1,
                         false,
                         this_table->pos,
                         ECS_TYPE_ecs_size_t,
                         comm);

  ecs_table__pos_en_regle(this_table);
}

/*----------------------------------------------------------------------------
 *  Fonction qui ecrit le contenu d'une table
 *   dans le fichier d'interface pour le Noyau
 *----------------------------------------------------------------------------*/

void
ecs_table_comm__ecr(ecs_table_t  *this_table,
                    const char   *comm_nom_rubrique,
                    size_t        location_id,
                    size_t        index_id,
                    size_t        n_location_values,
                    ecs_comm_t   *comm)
{
  size_t      nbr_val;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(this_table != NULL);

  nbr_val = ecs_table__ret_val_nbr(this_table);

  assert(nbr_val > 0);

  ecs_comm_write_section(comm_nom_rubrique,
                         nbr_val,
                         location_id,
                         index_id,
                         n_location_values,
                         false,
                         this_table->val,
                         ECS_TYPE_ecs_int_t,
                         comm);
}

/*----------------------------------------------------------------------------*/

