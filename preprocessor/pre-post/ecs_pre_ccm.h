#ifndef _ECS_PRE_CCM_H_
#define _ECS_PRE_CCM_H_

/*============================================================================
 * Read mesh from STAR-CCM+ format file
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

#include "cs_config.h"

#if defined(HAVE_CCM)

/*----------------------------------------------------------------------------
 * Fichiers `include' visibles du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/

#include "ecs_def.h"

/*----------------------------------------------------------------------------
 * Fichiers `include' visibles des paquetages visibles
 *----------------------------------------------------------------------------*/

#include "ecs_maillage.h"

/*----------------------------------------------------------------------------
 * Fichiers `include' publics  du  paquetage courant
 *----------------------------------------------------------------------------*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Read of mesh from a STAR-CCM+ format file.
 *
 * parameters:
 *   nom_fic_maillage <-- mesh file name
 *   num_maillage     <-- mesh number
 *----------------------------------------------------------------------------*/

ecs_maillage_t *
ecs_pre_ccm__lit_maillage(const char  *nom_fic_maillage,
                          int          num_maillage);

#endif /* HAVE_CCM */

/*----------------------------------------------------------------------------*/

#endif /* _ECS_PRE_CCM_H_ */

