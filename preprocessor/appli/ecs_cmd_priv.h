#ifndef _ECS_CMD_PRIV_H_
#define _ECS_CMD_PRIV_H_

/*============================================================================
 *  Définition privée de la structure `_ecs_cmd_t' decrivant
 *   les options de la ligne de commande
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


/*============================================================================
 *                                 Visibilité
 *============================================================================*/

#include "cs_config.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie standard C
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/

#include "ecs_def.h"
#include "ecs_tab_glob.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  des paquetages visibles
*----------------------------------------------------------------------------*/

#include "ecs_pre.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage courant
 *----------------------------------------------------------------------------*/


/*============================================================================
 *                       Définition de macros
 *============================================================================*/


#define ECS_CMD_LIGNE_MAX                                        1024


/*----------------------------------------------------------------------------
 * Définition des noms de fichiers et d'extensions
 *----------------------------------------------------------------------------*/

#define ECS_CMD_EXEC_NAME                             "cs_preprocess"

#define ECS_CMD_OUTFILE_NAME_DEFAULT               "preprocessor.log"

#define ECS_CMD_POST_CASE_DEFAULT                        "preprocess"
#define ECS_CMD_POST_DIR_ENS_EXT                           ".ensight"

/*----------------------------------------------------------------------------
 * Définition des mots-cles pour les options de la ligne de commande
 *----------------------------------------------------------------------------*/

#define ECS_CMD_KEY_MESH_GRP_SECTION                        "section"
#define ECS_CMD_KEY_MESH_GRP_ZONE                              "zone"

/*----------------------------------------------------------------------------
 *  Définition des options de la ligne de commande
 *----------------------------------------------------------------------------*/

#define ECS_CMD_OPTION_CASE                                  "--case"

#define ECS_CMD_OPTION_CWD                                    "--cwd"
#define ECS_CMD_OPTION_DUMP                                  "--dump"
#define ECS_CMD_OPTION_NULL_COMM                         "--no-write"
#define ECS_CMD_OPTION_FMT_MESH_FILE                       "--format"
#define ECS_CMD_OPTION_NUM_MESH                               "--num"
#define ECS_CMD_OPTION_GRP_CEL_MESH                       "--grp-cel"
#define ECS_CMD_OPTION_GRP_FAC_MESH                       "--grp-fac"

#define ECS_CMD_OPTION_HELP                                  "--help"
#define ECS_CMD_OPTION_HELP_1                                    "-h"

#if defined(HAVE_CGNS)
#define ECS_CMD_OPTION_POST_CGNS                             "--cgns"
#endif /* HAVE_CGNS */
#define ECS_CMD_OPTION_POST_ENS                           "--ensight"
#if defined(HAVE_MED)
#define ECS_CMD_OPTION_POST_MED                               "--med"
#endif /* HAVE_MED */

#define ECS_CMD_OPTION_POST_NO_POLY                  "--discard-poly"
#define ECS_CMD_OPTION_POST_TEXT                             "--text"
#define ECS_CMD_OPTION_POST_BIG_ENDIAN                 "--big-endian"
#define ECS_CMD_OPTION_POST_MAIN                           "--volume"
#define ECS_CMD_OPTION_POST_INFO                             "--info"

#define ECS_CMD_OPTION_MESH_FILE                             "--mesh"
#define ECS_CMD_OPTION_MESH_FILE_1                               "-m"
#define ECS_CMD_OPTION_ORIENT_CORREC                     "--reorient"
#define ECS_CMD_OPTION_OUTPUT_FILE                            "--log"

#define ECS_CMD_OPTION_VERSION                            "--version"

/*============================================================================
 *                       Définition des structures
 *============================================================================*/

struct _ecs_cmd_post_t {

  bool       no_poly;        /* suppression des polygones et polyedres */
  bool       simple;         /* pas de subdivision des "parts" (EnSight) */
  bool       text;           /* version texte (EnSight) */
  bool       big_endian;     /* binaire big-endian (EnSight) */

  bool       volume;         /* sortie du volume */
  bool       info;           /* sortie des maillages d'information */

};

struct _ecs_cmd_t {

  char                   *nom_cas;
  ecs_int_t               nbr_dump;
  ecs_int_t              *liste_num_maillage;
  ecs_pre_format_t       *liste_fmt_maillage;
  bool                   *liste_grp_maillage;
  ecs_tab_char_t          liste_fic_maillage;

#if defined(HAVE_CGNS)
  ecs_cmd_post_t         *post_cgns;
#endif /* HAVE_CGNS */
  ecs_cmd_post_t         *post_ens;
#if defined(HAVE_MED)
  ecs_cmd_post_t         *post_med;
#endif /* HAVE_MED */

  bool                    correct_orient;

  bool                    sim_comm;

};

/*----------------------------------------------------------------------------*/

#endif /* _ECS_CMD_PRIV_H_ */
