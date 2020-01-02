#ifndef _ECS_ELT_TYP_LISTE_H_
#define _ECS_ELT_TYP_LISTE_H_

/*============================================================================
 *  Définition d'un tableau constant statique
 *   définissant les types des éléments reconnus par le code
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "ecs_def.h"

BEGIN_C_DECLS

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct {

  int  elt_typ;
  int  som[4]; /* 4 sommets max. pour les sous-éléments */

} ecs_sous_elt_t;


typedef struct {

  char            nom[11]; /* Nom le plus long: polyhedron */
  int             nbr_som;
  int             nbr_sous_elt;
  ecs_sous_elt_t  sous_elt[6]; /* 6 sous-éléments max. */
} ecs_fic_elt_typ_t;

/*============================================================================
 * Global variable definitions
 *============================================================================*/

static const ecs_fic_elt_typ_t ecs_fic_elt_typ_liste_c[ECS_ELT_TYP_FIN] = {

  {                                                /* ECS_ELT_TYP_NUL       0 */
    "",
    0,
    0,
    {
      {0,{0}}
    }
  },

  {                                                /* ECS_ELT_TYP_FAC_TRIA  1 */
    "tria3",                                       /*                         */
    3,                                             /*        x 3              */
    0,                                             /*       / \               */
    {                                              /*      /   \              */
      {0, {0}}                                     /*     /     \             */
    }                                              /*  1 x-------x 2          */
  },

  {                                                /* ECS_ELT_TYP_FAC_QUAD  2 */
    "quad4",                                       /*                         */
    4,                                             /*  4 x-------x 3          */
    0,                                             /*    |       |            */
    {                                              /*    |       |            */
      {0, {0}}                                     /*    |       |            */
    }                                              /*  1 x-------x 2          */
  },

  {                                                /* ECS_ELT_TYP_CEL_TETRA 3 */
    "tetra4",                                      /*                         */
    4,                                             /*        x 4              */
    4,                                             /*       /|\               */
    {                                              /*      / | \              */
      {ECS_ELT_TYP_FAC_TRIA, { 1,  3,  2 }},       /*     /  |  \             */
      {ECS_ELT_TYP_FAC_TRIA, { 1,  2,  4 }},       /*  1 x- -|- -x 3          */
      {ECS_ELT_TYP_FAC_TRIA, { 1,  4,  3 }},       /*     \  |  /             */
      {ECS_ELT_TYP_FAC_TRIA, { 2,  3,  4 }}        /*      \ | /              */
    }                                              /*       \|/               */
  },                                               /*        x 2              */

  {                                                /* ECS_ELT_TYP_CEL_PYRAM 4 */
    "pyramid5",                                    /*                         */
    5,                                             /*         5 x             */
    5,                                             /*          /|\            */
    {                                              /*         //| \           */
      {ECS_ELT_TYP_FAC_TRIA, { 1,  2,  5 }    },   /*        // |  \          */
      {ECS_ELT_TYP_FAC_TRIA, { 1,  5,  4 }    },   /*     4 x/--|---x 3       */
      {ECS_ELT_TYP_FAC_TRIA, { 2,  3,  5 }    },   /*      //   |  /          */
      {ECS_ELT_TYP_FAC_TRIA, { 3,  4,  5 }    },   /*     //    | /           */
      {ECS_ELT_TYP_FAC_QUAD, { 1,  4,  3,  2 }}    /*  1 x-------x 2          */
    }
  },

  {                                                /* ECS_ELT_TYP_CEL_PRISM 5 */
    "penta6",                                      /*                         */
    6,                                             /*  4 x-------x 6          */
    5,                                             /*    |\     /|            */
    {                                              /*    | \   / |            */
      {ECS_ELT_TYP_FAC_TRIA, { 1,  3,  2 }    },   /*  1 x- \-/ -x 3          */
      {ECS_ELT_TYP_FAC_TRIA, { 4,  5,  6 }    },   /*     \ 5x  /             */
      {ECS_ELT_TYP_FAC_QUAD, { 1,  2,  5,  4 }},   /*      \ | /              */
      {ECS_ELT_TYP_FAC_QUAD, { 1,  4,  6,  3 }},   /*       \|/               */
      {ECS_ELT_TYP_FAC_QUAD, { 2,  3,  6,  5 }}    /*        x 2              */
    }
  },

  {                                                /* ECS_ELT_TYP_CEL_HEXA  6 */
    "hexa8",                                       /*                         */
    8,                                             /*                         */
    6,                                             /*     8 x-------x 7       */
    {                                              /*      /|      /|         */
      {ECS_ELT_TYP_FAC_QUAD, { 1,  4,  3,  2 }},   /*     / |     / |         */
      {ECS_ELT_TYP_FAC_QUAD, { 1,  2,  6,  5 }},   /*  5 x-------x6 |         */
      {ECS_ELT_TYP_FAC_QUAD, { 1,  5,  8,  4 }},   /*    | 4x----|--x 3       */
      {ECS_ELT_TYP_FAC_QUAD, { 2,  3,  7,  6 }},   /*    | /     | /          */
      {ECS_ELT_TYP_FAC_QUAD, { 3,  4,  8,  7 }},   /*    |/      |/           */
      {ECS_ELT_TYP_FAC_QUAD, { 5,  6,  7,  8 }}    /*  1 x-------x 2          */
    }
  },

  {                                                /* ECS_ELT_TYP_FAC_POLY  7 */
    "polygon",                                     /*                         */
    0,
    0,
    {
      {0,{0}}
    }
  },

  {                                                /* ECS_ELT_TYP_CEL_POLY  8 */
    "polyhedron",                                  /*                         */
    0,
    0,
    {
      {0,{0}}
    }
  }

};

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* _ECS_ELT_TYP_LISTE_H_ */
