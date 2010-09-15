#ifndef _ECS_CMD_PRIV_H_
#define _ECS_CMD_PRIV_H_

/*============================================================================
 * Definition of a `_ecs_cmd_t' structure tracking command-line options
 *============================================================================*/

/*
  This file is part of the Code_Saturne Preprocessor, element of the
  Code_Saturne CFD tool.

  Copyright (C) 1999-2010 EDF S.A., France

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

/*----------------------------------------------------------------------------*/

#include "cs_config.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "ecs_def.h"
#include "ecs_tab_glob.h"

#include "ecs_pre.h"

/*----------------------------------------------------------------------------
 * Headers for the current file
 *----------------------------------------------------------------------------*/

/*============================================================================
 * Type definitions
 *============================================================================*/

struct _ecs_cmd_post_t {

  bool       volume;         /* output volume mesh */
  bool       info;           /* output information meshes */

};

struct _ecs_cmd_t {

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

#if defined(HAVE_CGNS)
  ecs_cmd_post_t         *post_cgns;
#endif /* HAVE_CGNS */
  ecs_cmd_post_t         *post_ens;
#if defined(HAVE_MED)
  ecs_cmd_post_t         *post_med;
#endif /* HAVE_MED */

  bool                    correct_orient;

};

/*----------------------------------------------------------------------------*/

#endif /* _ECS_CMD_PRIV_H_ */
