#ifndef _ECS_CMD_H_
#define _ECS_CMD_H_

/*============================================================================
 *  Prototypes des fonctions de base
 *   associées à la structure `ecs_cmd_t' décrivant
 *   les options de la ligne de commande
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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

#include "cs_config.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "ecs_def.h"
#include "ecs_pre.h"

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct {

  char                   *fic_maillage;
  char                   *nom_cas;
  char                   *nom_out;
  int                     nbr_dump;
  int                     n_num_maillage;
  int                    *num_maillage;
  ecs_pre_format_t        fmt_maillage;
  bool                    grp_cel_section;
  bool                    grp_cel_zone;
  bool                    grp_fac_section;
  bool                    grp_fac_zone;

  bool                    correct_orient;

  char                    post_err[8];
  char                    post_vol[8];

} ecs_cmd_t;

/*============================================================================
 * Public function prototypes
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
