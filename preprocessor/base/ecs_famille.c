/*============================================================================
 *  Définitions des fonctions de base
 *   associées à la structure `ecs_famille_t' décrivant une famille
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


/*============================================================================
 *                                 Visibilité
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie standard C
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <string.h>  /* memcmp() */


/*----------------------------------------------------------------------------
 *  Fichiers `include' système
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/

#include "ecs_def.h"
#include "ecs_fic.h"
#include "ecs_file.h"
#include "ecs_mem.h"


/*----------------------------------------------------------------------------
  *  Fichiers `include' visibles des paquetages visibles
 *----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_descr_chaine.h"


/*----------------------------------------------------------------------------
 *  Fichier  `include' du  paquetage courant associe au fichier courant
 *----------------------------------------------------------------------------*/

#include "ecs_famille.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' prives   du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_famille_priv.h"


/*============================================================================
 *                              Fonctions privées
 *============================================================================*/


/*============================================================================
 *                             Fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *    Fonction de création d'une structure de famille `ecs_famille_t'
 *----------------------------------------------------------------------------*/

ecs_famille_t *
ecs_famille__cree(int           num,
                  ecs_descr_t  *descr_tete)
{
  ecs_famille_t * fam_loc;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Allocation et initialisation de la structure */
  /*==============================================*/

  /* Allocation de la structure globale de la famille */
  /*--------------------------------------------------*/

  ECS_MALLOC(fam_loc, 1, ecs_famille_t);

  /* Initialisation de la structure */
  /*--------------------------------*/

  fam_loc->num           = num;
  fam_loc->descr         = descr_tete;
  fam_loc->l_famille_sui = NULL;

  return fam_loc;
}

/*----------------------------------------------------------------------------
 *  Fonction libérant la structure `ecs_famille_t' donnée en argument.
 *  Elle renvoie un pointeur NULL
 *----------------------------------------------------------------------------*/

ecs_famille_t *
ecs_famille__detruit(ecs_famille_t  *this_fam)
{
  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(this_fam != NULL);

  ecs_descr_chaine__detruit(&this_fam->descr);

  ECS_FREE(this_fam);

  assert(this_fam == NULL);

  return this_fam;
}

/*----------------------------------------------------------------------------
 *  Fonction imprimant le contenu d'une structure `ecs_famille_t' donnée
 *   sur le flux décrit par la structure `ecs_file_t'
 *----------------------------------------------------------------------------*/

void
ecs_famille__imprime(const ecs_famille_t  *this_fam,
                     int                   imp_col,
                     FILE                 *fic_imp)
{
#define ECS_FCT_IMP_FAMILLE_NUM           "numero"
#define ECS_FCT_IMP_FAMILLE_DESCR         "descr"
#define ECS_FCT_IMP_FAMILLE_L_FAMILLE     "l_famille_sui"

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(this_fam != NULL);

  imp_col++;

  /* Impression du numero de famille */

  ecs_fic__imprime_val(fic_imp, imp_col, ECS_FCT_IMP_FAMILLE_NUM,
                       ECS_TYPE_ecs_int_t, &(this_fam->num));

  /* Impression de la liste chainee des descripteurs composant la famille */
  /*----------------------------------------------------------------------*/

  /* Impression du pointeur sur le descripteur de tete */

  ecs_fic__imprime_val(fic_imp, imp_col, ECS_FCT_IMP_FAMILLE_DESCR,
                       ECS_TYPE_void, this_fam->descr);

  /* Appel à la fonction d'impression d'une chaine de descripteurs */

  ecs_descr_chaine__imprime(this_fam->descr,
                            imp_col + 1,
                            fic_imp);

  /* Impression du lien sur une eventuelle famille suivante */
  /*--------------------------------------------------------*/

  ecs_fic__imprime_val(fic_imp, imp_col, ECS_FCT_IMP_FAMILLE_L_FAMILLE,
                       ECS_TYPE_void, this_fam->l_famille_sui);
}

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie la taille en octets d'une structure `ecs_famille_t'
 *----------------------------------------------------------------------------*/

float
ecs_famille__ret_taille(const ecs_famille_t  *this_fam)
{
  float  taille;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(this_fam != NULL);

  taille = sizeof(*this_fam);

  taille += ecs_descr_chaine__ret_taille(this_fam->descr);

  return taille;
}

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie le nombre et la liste des des pointeurs sur les noms
 *   des descripteurs de la famille donnée en argument
 *----------------------------------------------------------------------------*/

ecs_tab_char_t
ecs_famille__ret_nom(const ecs_famille_t  *this_fam)
{
  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(this_fam != NULL);

  return ecs_descr_chaine__ret_nom(this_fam->descr);
}

/*----------------------------------------------------------------------------
 *  Fonction qui alloue une structure `ecs_famille_t' et qui remplit
 *   son contenu en copiant le contenu de la structure donnée en argument
 *----------------------------------------------------------------------------*/

ecs_famille_t *
ecs_famille__copie(ecs_famille_t  *this_famille)
{
  ecs_descr_t     * descr;
  ecs_famille_t   * famille_loc;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(this_famille != NULL);

  descr = ecs_descr_chaine__copie(this_famille->descr);

  famille_loc = ecs_famille__cree(this_famille->num, descr);

  return famille_loc;
}

/*----------------------------------------------------------------------------
 *  Fonction qui affiche la définition de la famille
 *----------------------------------------------------------------------------*/

void
ecs_famille__affiche(const ecs_famille_t  *this_fam)
{
  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(this_fam != NULL);

  printf("  %s %d\n",
         _("Family"), this_fam->num);

  ecs_descr_chaine__affiche(this_fam->descr,
                            (int)strlen(_("Family")) + 1);
}

/*----------------------------------------------------------------------------*/

