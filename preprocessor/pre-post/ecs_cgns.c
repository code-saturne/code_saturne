/*============================================================================
 * Définitions de base pour le support CGNS
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


#include "cs_config.h"

#if defined(HAVE_CGNS)


/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie CGNS
 *----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif

#include <cgnslib.h>

#ifdef __cplusplus
}
#endif


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/

#include "ecs_def.h"


/*----------------------------------------------------------------------------
 *  Fichier  `include' du  paquetage courant associé au fichier courant
 *----------------------------------------------------------------------------*/

#include "ecs_cgns.h"


/*============================================================================
 * Définitions de variables globales
 *============================================================================*/

/* Définition des éléments */
/*=========================*/

/* Le tableau suivant est donné dans le même ordre que les définitions */
/* de cgnslib.h (mais ne contient pas les deux premières entrées)      */

const ecs_cgns_elt_t
ecs_cgns_elt_liste_c[ECS_CGNS_NBR_TYP_ELT] = {
  {
                               /* 1 */
    CS_CG_ENUM(NODE),
    ECS_ELT_TYP_NUL,
    1,
    { 0 }
  },
  {                                /* 2 */
    CS_CG_ENUM(BAR_2),
    ECS_ELT_TYP_NUL,
    2,
    { 0 }
  },
  {                                /* 3 */
    CS_CG_ENUM(BAR_3),
    ECS_ELT_TYP_NUL,
    3,
    { 0 }
  },
  {                                /* 4 */
    CS_CG_ENUM(TRI_3),
    ECS_ELT_TYP_FAC_TRIA,
    3,
    { 1, 2, 3 }
  },
  {                                /* 5 */
    CS_CG_ENUM(TRI_6),
    ECS_ELT_TYP_FAC_TRIA,
    6,
    { 1, 2, 3 }
  },
  {                                /* 6 */
    CS_CG_ENUM(QUAD_4),
    ECS_ELT_TYP_FAC_QUAD,
    4,
    { 1, 2, 3, 4 }
  },
  {                                /* 7 */
    CS_CG_ENUM(QUAD_8),
    ECS_ELT_TYP_FAC_QUAD,
    8,
    { 1, 2, 3, 4 }
  },
  {                                /* 8 */
    CS_CG_ENUM(QUAD_9),
    ECS_ELT_TYP_FAC_QUAD,
    9,
    { 1, 2, 3, 4 }
  },
  {                                /* 9 */
    CS_CG_ENUM(TETRA_4),
    ECS_ELT_TYP_CEL_TETRA,
    4,
    { 1, 2, 3, 4 }
  },
  {                               /* 10 */
    CS_CG_ENUM(TETRA_10),
    ECS_ELT_TYP_CEL_TETRA,
    10,
    { 1, 2, 3, 4 }
  },
  {                               /* 11 */
    CS_CG_ENUM(PYRA_5),
    ECS_ELT_TYP_CEL_PYRAM,
    5,
    { 1, 2, 3, 4, 5 }
  },
  {                               /* 12 */
    CS_CG_ENUM(PYRA_14),
    ECS_ELT_TYP_CEL_PYRAM,
    14,
    { 1, 2, 3, 4, 5 }
  },
  {                               /* 13 */
    CS_CG_ENUM(PENTA_6),
    ECS_ELT_TYP_CEL_PRISM,
    6,
    { 1, 2, 3, 4, 5, 6 }
  },
  {                               /* 14 */
    CS_CG_ENUM(PENTA_15),
    ECS_ELT_TYP_CEL_PRISM,
    15,
    { 1, 2, 3, 4, 5, 6 }
  },
  {                               /* 15 */
    CS_CG_ENUM(PENTA_18),
    ECS_ELT_TYP_CEL_PRISM,
    18,
    { 1, 2, 3, 4, 5, 6 }
  },
  {                               /* 16 */
    CS_CG_ENUM(HEXA_8),
    ECS_ELT_TYP_CEL_HEXA,
    8,
    { 1, 2, 3, 4, 5, 6, 7, 8 },
  },
  {                               /* 17 */
    CS_CG_ENUM(HEXA_20),
    ECS_ELT_TYP_CEL_HEXA,
    20,
    { 1, 2, 3, 4, 5, 6, 7, 8 },
  },
  {                               /* 18 */
    CS_CG_ENUM(HEXA_27),
    ECS_ELT_TYP_CEL_HEXA,
    27,
    { 1, 2, 3, 4, 5, 6, 7, 8 },
  },
  {                               /* 19 */
    CS_CG_ENUM(MIXED),
    ECS_ELT_TYP_NUL,
    0,
    { 0 },
  },
  {                               /* 20 */
    CS_CG_ENUM(NGON_n),
    ECS_ELT_TYP_FAC_POLY,
    0,
    { 0 }
  }
};

#endif /* HAVE_CGNS */

/*----------------------------------------------------------------------------*/

