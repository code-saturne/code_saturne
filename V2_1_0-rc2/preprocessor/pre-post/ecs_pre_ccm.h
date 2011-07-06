/*============================================================================
 *
 *     This file is part of the Code_Saturne Preprocessor, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2009 EDF S.A., France
 *
 *     contact: saturne-support@edf.fr
 *
 *     The Code_Saturne Kernel is free software; you can redistribute it
 *     and/or modify it under the terms of the GNU General Public License
 *     as published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version.
 *
 *     The Code_Saturne Kernel is distributed in the hope that it will be
 *     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with the Code_Saturne Kernel; if not, write to the
 *     Free Software Foundation, Inc.,
 *     51 Franklin St, Fifth Floor,
 *     Boston, MA  02110-1301  USA
 *
 *============================================================================*/

#ifndef _ECS_PRE_CCM_H_
#define _ECS_PRE_CCM_H_

/*============================================================================
 * Read mesh from STAR-CCM+ format file
 *============================================================================*/

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

