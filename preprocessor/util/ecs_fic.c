/*============================================================================
 *  Définitions des fonctions de base
 *   associées aux impressions dans un fichier
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

#include "ecs_def.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdarg.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "ecs_file.h"
#include "ecs_mem.h"

/*----------------------------------------------------------------------------
 *  Headers for the current file
 *----------------------------------------------------------------------------*/

#include "ecs_fic.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/* Formats pour impressions de tableaux dans un fichier */

#define ECS_LOC_IMP_FMT_INT             "%12d"
#define ECS_LOC_IMP_FMT_REAL         "%#12.5E"
#define ECS_LOC_IMP_FMT_CHAR            "%12s"
#define ECS_LOC_IMP_FMT_SIZE_T         "%12lu"
#define ECS_LOC_IMP_FMT_PTR             "%12p"
#define ECS_LOC_IMP_FMT_CONT    "............"
#define ECS_LOC_IMP_FMT_WIDTH               12

#define ECS_LOC_IMP_INDENT                   8

#define ECS_LOC_IMP_COL_VAL                  9
#define ECS_LOC_IMP_COL_IND                 51

/*============================================================================
 * Local function definitions
 *============================================================================*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction d'impression du nom et de la valeur d'un pointeur
 *----------------------------------------------------------------------------*/

void
ecs_fic__imprime_ptr(FILE         *fic_imp,
                     const int     profondeur_imp,
                     const char   *nom_ptr,
                     const void   *ptr)
{
  /* Instructions */

  assert(fic_imp != NULL);

  fprintf(fic_imp,
          "%-*s%s  & ", ECS_LOC_IMP_INDENT * profondeur_imp, "",
          nom_ptr);

  if (ptr == NULL)
    fprintf(fic_imp, "%s\n", "NULL");

  else
    fprintf(fic_imp, "%p\n", ptr);
}

/*----------------------------------------------------------------------------
 *  Fonction d'impression du nom et de la valeur d'une variable
 *----------------------------------------------------------------------------*/

void
ecs_fic__imprime_val(FILE         *fic_imp,
                     const int     profondeur_imp_nom,
                     const char   *nom,
                     ecs_type_t    typ_e,
                     const void   *val)
{
  /* Variables locales */

  int  nb_col_nom;
  int  profondeur_imp_val;

  /* Instructions */

  assert(fic_imp != NULL);

  if ((int)strlen(nom) < ECS_LOC_IMP_INDENT)
    nb_col_nom = 1;
  else
    nb_col_nom = 2;

  profondeur_imp_val
    = ECS_LOC_IMP_COL_VAL - profondeur_imp_nom - nb_col_nom - 1;

  fprintf(fic_imp,
          "%-*s%-*s%-*s",
          ECS_LOC_IMP_INDENT * profondeur_imp_nom, "" ,
          ECS_LOC_IMP_INDENT * nb_col_nom        , nom,
          ECS_LOC_IMP_INDENT * profondeur_imp_val, "");

  if (val == NULL)

    fprintf(fic_imp, ECS_LOC_IMP_FMT_CHAR "\n", "NULL");

  else {

    switch(typ_e) {

    case ECS_TYPE_ecs_int_t:

      {
        const int val_imp = *(const int *)val;
        fprintf(fic_imp, ECS_LOC_IMP_FMT_INT "\n", val_imp);
      }
      break;

    case ECS_TYPE_ecs_coord_t:

      {
        const ecs_coord_t val_imp = *(const ecs_coord_t *)val;
        fprintf(fic_imp, ECS_LOC_IMP_FMT_REAL "\n", val_imp);
      }
      break;

    case ECS_TYPE_char:

      {
        const char *val_imp = (const char *)val;
        fprintf(fic_imp, ECS_LOC_IMP_FMT_CHAR "\n", val_imp);
      }
      break;

    case ECS_TYPE_size_t:

      {
        const size_t val_imp = *(const size_t *)val;
        fprintf(fic_imp, ECS_LOC_IMP_FMT_SIZE_T "\n", (unsigned long)val_imp);
      }
      break;

    case ECS_TYPE_void:

      fprintf(fic_imp, ECS_LOC_IMP_FMT_PTR  "\n", val);
      break;

    default:
      assert(   typ_e == ECS_TYPE_ecs_int_t
             || typ_e == ECS_TYPE_ecs_coord_t
             || typ_e == ECS_TYPE_char
             || typ_e == ECS_TYPE_size_t
             || typ_e == ECS_TYPE_void);
    }

  }

  fflush(fic_imp);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
