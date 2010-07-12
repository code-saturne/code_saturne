#ifndef _ECS_CHAMP_PRIV_H_
#define _ECS_CHAMP_PRIV_H_

/*============================================================================
 *  Définition privée de la structure `_ecs_champ_t' décrivant un champ
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

#include "ecs_champ.h"


/*============================================================================
 *                       Définition de macros
 *============================================================================*/


/*============================================================================
 *                       Définition des structures
 *============================================================================*/

struct _ecs_champ_t {

  size_t                nbr;           /* Nombre d'éléments associés */

  size_t                pas;           /* Pas des positions si elles
                                          forment une REGLE (0 sinon) */
  ecs_size_t           *pos;           /* Positions des valeurs du champ
                                          (NULL si elles forment une REGLE) */

  ecs_int_t            *val;           /* Valeurs du champ */

  ecs_descr_t          *descr;         /* Tête de la liste chaînée des
                                          descripteurs de champ "attribut" */
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

#endif /* _ECS_CHAMP_PRIV_H_ */
