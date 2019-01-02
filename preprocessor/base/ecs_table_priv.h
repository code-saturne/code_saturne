#ifndef _ECS_CHAMP_PRIV_H_
#define _ECS_CHAMP_PRIV_H_

/*============================================================================
 *  Définition privée de la structure `_ecs_table_t' décrivant une table
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie standard C
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  des paquetages visibles
 *----------------------------------------------------------------------------*/

#include "ecs_descr.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_table.h"


/*============================================================================
 *                       Définition de macros
 *============================================================================*/


/*============================================================================
 *                       Définition des structures
 *============================================================================*/

struct _ecs_table_t {

  size_t                nbr;           /* Nombre d'éléments associés */

  size_t                pas;           /* Pas des positions si elles
                                          forment une REGLE (0 sinon) */
  ecs_size_t           *pos;           /* Positions des valeurs de la table
                                          (NULL si elles forment une REGLE) */

  ecs_int_t            *val;           /* Valeurs de la table */

  ecs_descr_t          *descr;         /* Tête de la liste chaînée des
                                          descripteurs de table "attribut" */
};

/*============================================================================
 *
 *  Schéma d'association entre les tables "positions" et "valeurs"
 * ----------------------------------------------------------------
 *
 * Table des valeurs `val' (de dimension `pos[nbr-1]-1')
 *
 * .---.-------..------.-------.------..------.-------..-------.------.
 * |   |  ...  ||      |  ...  |      ||      |  ...  ||  ...  |      |
 * `---'-------'`------'-------'------'`------'-------'`-------'------'
 *   0           iVal-1         jVal-2  jVal-1                  nVal-2 nVal-1
 *
 *                  |                      |                              |
 *                  |                      |                              |
 *                  `----------.       .---'          .-------------------'
 *                             |       |              |
 *            .-----.-------.------.------.-------.------.
 *            |  1  |  ...  | iVal | jVal |  ...  | nVal |
 *            `-----'-------'------'------'-------'------'
 *               0            iPos  iPos+1          nPos = pos->nbr - 1
 *
 * Table des positions `pos' (de dimension `nbr' + 1)
 *
 *============================================================================*/

/*----------------------------------------------------------------------------*/

#endif /* _ECS_TABLE_PRIV_H_ */
