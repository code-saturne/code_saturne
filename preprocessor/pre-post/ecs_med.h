#ifndef _ECS_MED_H_
#define _ECS_MED_H_

/*============================================================================
 *  Définitions de macros et constantes
 *   servant au traitement de maillages au format MED
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

#if defined(HAVE_MED)


/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie standard C
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/

#include "ecs_def.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  des paquetages visibles
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage courant
 *----------------------------------------------------------------------------*/


/*============================================================================
 *                      Définitions de variables globales
 *============================================================================*/

/* Variables de définition de la version des librairies */

extern int  ecs_glob_med_ver_maj;
extern int  ecs_glob_med_ver_min;
extern int  ecs_glob_med_ver_rel;
extern int  ecs_glob_hdf5_ver_maj;
extern int  ecs_glob_hdf5_ver_min;
extern int  ecs_glob_hdf5_ver_rel;


/*============================================================================
 *                          Définitions de types
 *============================================================================*/

/* Structure définissant un cas MED */

typedef struct _ecs_med_t ecs_med_t;


/*============================================================================
 *                       Définition de macros
 *============================================================================*/


/*============================================================================
 * Fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Fonction de récupération d'informations sur les librairies associées
 * (pour les librairies chargées dynamiquement dont la version peut être
 * modifiée au lancement).
 *----------------------------------------------------------------------------*/

void
ecs_med__version_shlib(void);

#endif /* HAVE_MED */

/*----------------------------------------------------------------------------*/

#endif /* _ECS_MED_H_ */
