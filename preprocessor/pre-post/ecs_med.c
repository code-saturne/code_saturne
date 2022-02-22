/*============================================================================
 * Définitions de base pour le support MED
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


/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie standard C
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "cs_config.h"

/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie MED
 *----------------------------------------------------------------------------*/

#if defined(HAVE_MED)

#include "ecs_med_priv.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/

#include "ecs_def.h"


/*----------------------------------------------------------------------------
 *  Fichier  `include' du  paquetage courant associé au fichier courant
 *----------------------------------------------------------------------------*/

#include "ecs_med.h"
#include "ecs_med_priv.h"


/*=============================================================================
 * Définitions de variables globales
 *============================================================================*/

/* Support de la librairie d'échange de maillages et champs MED */

#if defined(MED_NUM_MAJEUR)
int  ecs_glob_med_ver_maj  = MED_NUM_MAJEUR;
int  ecs_glob_med_ver_min  = MED_NUM_MINEUR;
int  ecs_glob_med_ver_rel  = MED_NUM_RELEASE;
#else
int  ecs_glob_med_ver_maj  = 2;
int  ecs_glob_med_ver_min  = 3;
int  ecs_glob_med_ver_rel  = -1;
#endif
int  ecs_glob_hdf5_ver_maj = H5_VERS_MAJOR;
int  ecs_glob_hdf5_ver_min = H5_VERS_MINOR;
int  ecs_glob_hdf5_ver_rel = H5_VERS_RELEASE;


const ecs_fic_med_init_elt_t
ecs_fic_med_init_elt_liste_c[ECS_MED_NBR_TYP_ELT] = {
  {                                /* 1 */
    MED_TRIA3,
    ECS_ELT_TYP_FAC_TRIA,
    ECS_MED_ORDER_LINEAR,
    { 1, 2, 3 }
  },
  {                                /* 2 */
    MED_TRIA6,
    ECS_ELT_TYP_FAC_TRIA,
    ECS_MED_ORDER_PARABOLIC,
    { 1, 2, 3 }
  },
  {                                /* 3 */
    MED_QUAD4,
    ECS_ELT_TYP_FAC_QUAD,
    ECS_MED_ORDER_LINEAR,
    { 1, 2, 3, 4 }
  },
  {                                /* 4 */
    MED_QUAD8,
    ECS_ELT_TYP_FAC_QUAD,
    ECS_MED_ORDER_PARABOLIC,
    { 1, 2, 3, 4 }
  },
  {                                /* 5 */
    MED_TETRA4,
    ECS_ELT_TYP_CEL_TETRA,
    ECS_MED_ORDER_LINEAR,
    { 1, 3, 2, 4 }
  },
  {                                /* 6 */
    MED_TETRA10,
    ECS_ELT_TYP_CEL_TETRA,
    ECS_MED_ORDER_PARABOLIC,
    { 1, 3, 2, 4 }
  },
  {                               /*  7 */
    MED_PYRA5,
    ECS_ELT_TYP_CEL_PYRAM,
    ECS_MED_ORDER_LINEAR,
    { 1, 4, 3, 2, 5 }
  },
  {                               /*  8 */
    MED_PYRA13,
    ECS_ELT_TYP_CEL_PYRAM,
    ECS_MED_ORDER_PARABOLIC,
    { 1, 4, 3, 2, 5 }
  },
  {                               /*  9 */
    MED_PENTA6,
    ECS_ELT_TYP_CEL_PRISM,
    ECS_MED_ORDER_LINEAR,
    { 1, 3, 2, 4, 6, 5 }
  },
  {                               /* 10 */
    MED_PENTA15,
    ECS_ELT_TYP_CEL_PRISM,
    ECS_MED_ORDER_PARABOLIC,
    { 1, 3, 2, 4, 6, 5 }
  },
  {                               /* 11 */
    MED_HEXA8,
    ECS_ELT_TYP_CEL_HEXA,
    ECS_MED_ORDER_LINEAR,
    { 1, 4, 3, 2, 5, 8, 7, 6 },
  },
  {                               /* 12 */
    MED_HEXA20,
    ECS_ELT_TYP_CEL_HEXA,
    ECS_MED_ORDER_PARABOLIC,
    { 1, 4, 3, 2, 5, 8, 7, 6 },
  },
  {                               /* 13 */
    MED_POLYGON,
    ECS_ELT_TYP_FAC_POLY,
    ECS_MED_ORDER_LINEAR,
    { 0 },
  },
  {                               /* 14 */
    MED_POLYHEDRON,
    ECS_ELT_TYP_CEL_POLY,
    ECS_MED_ORDER_LINEAR,
    { 0 }
  }
};


/*============================================================================
 * Fonctions privées
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
ecs_med__version_shlib(void)
{
  med_int   med_majeur;
  med_int   med_mineur;
  med_int   med_release;

  MEDlibraryNumVersion(&med_majeur, &med_mineur, &med_release);

  ecs_glob_med_ver_maj  = med_majeur;
  ecs_glob_med_ver_min  = med_mineur;
  ecs_glob_med_ver_rel  = med_release;

  {
    med_int  hdf5_majeur;
    med_int  hdf5_mineur;
    med_int  hdf5_release;

    MEDlibraryHdfNumVersion(&hdf5_majeur, &hdf5_mineur, &hdf5_release);

    ecs_glob_hdf5_ver_maj = hdf5_majeur;
    ecs_glob_hdf5_ver_min = hdf5_mineur;
    ecs_glob_hdf5_ver_rel = hdf5_release;
  }
}

/*----------------------------------------------------------------------------*/

#endif /* HAVE_MED */
