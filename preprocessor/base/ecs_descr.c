/*============================================================================
 * Definitions of base functions associated to the `ecs_descr_t' structure
 * describing an array (ecs_table_t) descriptor.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "ecs_def.h"
#include "ecs_fic.h"
#include "ecs_mem.h"

/*----------------------------------------------------------------------------
 *  Headers for the current file
 *----------------------------------------------------------------------------*/

#include "ecs_descr.h"
#include "ecs_descr_priv.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local function defintions
 *============================================================================*/

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*=============================================================================
 * Local type definitions
 *============================================================================*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction de création d'une structure de descripteur de table
 *----------------------------------------------------------------------------*/

ecs_descr_t *
ecs_descr__cree(int          ide,
                const char  *nom)
{
  ecs_descr_t * descr_loc;

  /* Allocation de la structure globale du descripteur de table */

  ECS_MALLOC(descr_loc, 1, ecs_descr_t);

  /* Initialisation du numéro du descripteur */

  descr_loc->num = 1;

  /* Affectation du nom du descripteur */

  if (nom != NULL) {
    ECS_MALLOC(descr_loc->nom, strlen(nom) + 1, char);
    strcpy(descr_loc->nom, nom);
  }
  else {
    int _tmp_ide = ECS_ABS(ide);
    int l = 1;
    while (_tmp_ide/10 > 0) {
      _tmp_ide /= 10;
      l += 1;
    }
    if (ide < 0)
      l += 1;
    ECS_MALLOC(descr_loc->nom, l + 1, char);
    sprintf(descr_loc->nom, "%d", ide);
  }

  /* Initialisation par défaut du lien sur le descripteur suivant */

  descr_loc->l_descr_sui = NULL;

  /* On renvoie un pointeur sur la structure */

  return descr_loc;
}

/*----------------------------------------------------------------------------
 *  Fonction libérant la structure `ecs_descr_t' donnée en argument.
 *  Elle renvoie un pointeur NULL
 *----------------------------------------------------------------------------*/

ecs_descr_t *
ecs_descr__detruit(ecs_descr_t  *this_descr)
{
  assert(this_descr != NULL);

  /* Libération du nom du descripteur */

  if (this_descr->nom != NULL)
    ECS_FREE(this_descr->nom);

  /* Libération de la structure */

  ECS_FREE(this_descr);

  return this_descr;
}

/*----------------------------------------------------------------------------
 *  Fonction imprimant le contenu d'une structure `ecs_descr_t' donnée
 *   sur le flux décrit par la structure `FILE'
 *----------------------------------------------------------------------------*/

void
ecs_descr__imprime(const ecs_descr_t  *this_descr,
                   ecs_int_t           imp_col,
                   FILE               *fic_imp)
{
#define ECS_FCT_IMP_DESCR_NUM         "numero"
#define ECS_FCT_IMP_DESCR_NOM         "nom"
#define ECS_FCT_IMP_DESCR_L_DESCR     "l_descr_sui"

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(this_descr != NULL);

  imp_col++;

  /* Écriture du contenu de la structure `ecs_descr_t' */
  /*===============================================*/

  ecs_fic__imprime_val(fic_imp, imp_col, ECS_FCT_IMP_DESCR_NUM,
                       ECS_TYPE_ecs_int_t, &this_descr->num);

  ecs_fic__imprime_val(fic_imp, imp_col, ECS_FCT_IMP_DESCR_NOM,
                       ECS_TYPE_char, this_descr->nom);

  ecs_fic__imprime_val(fic_imp, imp_col, ECS_FCT_IMP_DESCR_L_DESCR,
                       ECS_TYPE_void, this_descr->l_descr_sui);
}

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie la taille en octets d'une structure `ecs_descr_t'
 *----------------------------------------------------------------------------*/

float
ecs_descr__ret_taille(const ecs_descr_t  *this_descr)
{
  float   taille;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(this_descr != NULL);

  taille = (float)sizeof(*this_descr);
  taille += (float)sizeof(*(this_descr->nom));

  return taille;
}

/*----------------------------------------------------------------------------
 *  Fonction qui alloue une structure `ecs_descr_t' et qui remplit
 *   son contenu en copiant le contenu de la structure donnée en argument
 *   sauf le numéro de descripteur qui n'est pas recopie
 *----------------------------------------------------------------------------*/

ecs_descr_t *
ecs_descr__copie(ecs_descr_t  *this_descr)
{
  ecs_descr_t   * descr_loc;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(this_descr != NULL);

  descr_loc = ecs_descr__cree(ECS_DESCR_IDE_NUL,
                              this_descr->nom);

  return descr_loc;
}

/*----------------------------------------------------------------------------
 *  Fonction qui compare 2 descripteurs
 *
 *  La fonction renvoie :
 *  - `true'  si les deux descripteurs ont les mêmes identificateur et nom
 *  - `false' sinon
 *----------------------------------------------------------------------------*/

bool
ecs_descr__compare(const ecs_descr_t  *descr_1,
                   const ecs_descr_t  *descr_2)
{
  bool       bool_retour = true;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(descr_1 != NULL && descr_2 != NULL);

  if (descr_1->nom != NULL || descr_2->nom != NULL) {

    if((descr_1->nom == NULL && descr_2->nom != NULL) ||
       (descr_1->nom != NULL && descr_2->nom == NULL) ||
       strcmp(descr_1->nom, descr_2->nom ) != 0 ) {

      bool_retour = false;

    }

  }

  return bool_retour;
}

/*----------------------------------------------------------------------------
 *  Fonction qui affiche le nom d'un descripteur
 *----------------------------------------------------------------------------*/

void
ecs_descr__affiche(const ecs_descr_t  *descr,
                   int                 decal)
{
  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(descr->nom != NULL);

  printf("  %*s%s \"%s\"\n",
         decal, "", _("Group"),
         descr->nom);
}

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie le nom du descripteur donné en argument
 *----------------------------------------------------------------------------*/

const char *
ecs_descr__ret_nom(const ecs_descr_t  *descr)
{
  assert(descr != NULL);

  return descr->nom;
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

