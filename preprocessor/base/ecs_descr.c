/*============================================================================
 *  Définitions des fonctions de base
 *   associées à la structure `ecs_descr_t' décrivant un descripteur de champ
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


/*============================================================================
 *                                 Visibilité
 *============================================================================*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie standard C
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <string.h>  /* memcmp() */


/*----------------------------------------------------------------------------
 *  Fichiers `include' système
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/

#include "ecs_def.h"
#include "ecs_fic.h"
#include "ecs_mem.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles des paquetages visibles
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage courant
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichier  `include' du  paquetage courant associé au fichier courant
 *----------------------------------------------------------------------------*/

#include "ecs_descr.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' privés   du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_descr_priv.h"


/*============================================================================
 *                              Fonctions privées
 *============================================================================*/


/*============================================================================
 *                             Fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *    Fonction de création d'une structure de descripteur de champ
 *----------------------------------------------------------------------------*/

ecs_descr_t *
ecs_descr__cree(ecs_descr_typ_t   typ,
                const ecs_int_t   ide,
                const char       *nom)
{
  ecs_descr_t * descr_loc;

  /* Allocation de la structure globale du descripteur de champ */

  ECS_MALLOC(descr_loc, 1, ecs_descr_t);

  /* Initialisation du numéro du descripteur */

  descr_loc->num = 1;

  /* Affectation du type du descripteur */

  descr_loc->typ = typ;

  /* Affectation de l'identificateur du descripteur */

  descr_loc->ide = ide;

  /* Affectation du nom du descripteur */

  if (nom != NULL) {
    ECS_MALLOC(descr_loc->nom, strlen(nom) + 1, char);
    strcpy(descr_loc->nom, nom);
  }
  else {
    descr_loc->nom = NULL;
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
#define ECS_FCT_IMP_DESCR_TYP         "type"
#define ECS_FCT_IMP_DESCR_IDE         "ide"
#define ECS_FCT_IMP_DESCR_NOM         "nom"
#define ECS_FCT_IMP_DESCR_L_DESCR     "l_descr_sui"

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(this_descr != NULL);

  imp_col++;

  /* Écriture du contenu de la structure `ecs_descr_t' */
  /*===============================================*/

  /* Numéro du descripteur */
  /*-----------------------*/

  ecs_fic__imprime_val(fic_imp, imp_col, ECS_FCT_IMP_DESCR_NUM,
                       ECS_TYPE_ecs_int_t, &this_descr->num);

  /* Type du descripteur */
  /*---------------------*/

  ecs_fic__imprime_val(fic_imp, imp_col, ECS_FCT_IMP_DESCR_TYP,
                       ECS_TYPE_ecs_int_t, &this_descr->typ);

  /* Identificateur du descripteur */
  /*-------------------------------*/

  ecs_fic__imprime_val(fic_imp, imp_col, ECS_FCT_IMP_DESCR_IDE,
                       ECS_TYPE_ecs_int_t, &this_descr->ide);

  /* Nom du descripteur */
  /*--------------------*/

  ecs_fic__imprime_val(fic_imp, imp_col, ECS_FCT_IMP_DESCR_NOM,
                       ECS_TYPE_char, this_descr->nom);

  /* Lien sur un éventuel descripteur */
  /*----------------------------------*/

  ecs_fic__imprime_val(fic_imp, imp_col, ECS_FCT_IMP_DESCR_L_DESCR,
                       ECS_TYPE_void, this_descr->l_descr_sui);
}

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie la taille en octets d'une structure `ecs_descr_t'
 *----------------------------------------------------------------------------*/

float
ecs_descr__ret_taille(const ecs_descr_t *this_descr)
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
ecs_descr__copie(ecs_descr_t * this_descr)
{
  ecs_descr_t   * descr_loc;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(this_descr != NULL);

  descr_loc = ecs_descr__cree(this_descr->typ,
                              this_descr->ide,
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
  bool       bool_retour;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(descr_1 != NULL && descr_2 != NULL);


  if (descr_1->typ == descr_2->typ && descr_1->ide == descr_2->ide)
    bool_retour = true;
  else
    bool_retour = false;


  if (bool_retour == true) {

    if (descr_1->nom != NULL || descr_2->nom != NULL) {

      if((descr_1->nom == NULL && descr_2->nom != NULL) ||
         (descr_1->nom != NULL && descr_2->nom == NULL) ||
         strcmp(descr_1->nom,
                descr_2->nom ) != 0 ) {

        bool_retour = false;

      }

    }

  }

  return bool_retour;
}

/*----------------------------------------------------------------------------
 *  Fonction qui compare 2 descripteurs
 *
 *  La fonction renvoie :
 *  - `true' si les types des deux identificateurs sont identiques et
 *           -       si les noms            des 2 descripteurs sont a `NULL'
 *                et si les identificateurs des 2 descripteurs sont identiques
 *           - ou    si les noms            des 2 descripteurs sont identiques
 *                   (et tous deux differents de `NULL')
 *  - `false' sinon
 *----------------------------------------------------------------------------*/

bool
ecs_descr__compare_selection(const ecs_descr_t  *descr_1,
                             const ecs_descr_t  *descr_2)
{
  assert(descr_1 != NULL && descr_2 != NULL);

  if (descr_1->typ != descr_2->typ)
    return false;

  if (descr_1->nom == NULL && descr_2->nom == NULL) {

    /* Les noms des 2 descripteurs sont a `NULL' : */
    /* on compare leurs identificateurs            */

    if (descr_1->ide != descr_2->ide)

      return false;

  }
  else {

    /* Au moins un nom des deux descripteurs n'est pas a `NULL' : */
    /*    on compare les deux noms                                */


    if((descr_1->nom == NULL && descr_2->nom != NULL) ||
       (descr_1->nom != NULL && descr_2->nom == NULL) ||
       strcmp(descr_1->nom,
              descr_2->nom ) != 0 )

      return false;
  }

  return true;
}

/*----------------------------------------------------------------------------
 *  Fonction qui affiche le contenu d'un descripteur
 *
 *  Sont affichés :
 *   - le type du descripteur ("Couleur" ou "Groupe")
 *   - l'identificateur dans le cas d'une couleur
 *   - la chaîne de caractères dans le cas d'un groupe
 *----------------------------------------------------------------------------*/

void
ecs_descr__affiche(const ecs_descr_t  *descr,
                   int                 decal)
{
  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  switch (descr->typ) {

  case ECS_DESCR_COULEUR:

    if (descr->nom == NULL)
      printf("  %*s%s %d\n",
             decal, "", _("Color"),
             descr->ide);
    else
      printf("  %*s%s %d (%s)\n",
             decal, "", _("Color"), descr->ide, descr->nom);

    break;

  case ECS_DESCR_GROUPE:

    assert(descr->nom != NULL);

    printf("  %*s%s \"%s\"\n",
           decal, "", _("Group"),
           descr->nom);

    break;

  default:

    assert(descr->typ == ECS_DESCR_COULEUR ||
           descr->typ == ECS_DESCR_GROUPE     );

  }
}

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie le type du descripteur :
 *  - soit `ECS_DESCR_TYP_COULEUR' pour un descripteur de type "couleur"
 *  - soit `ECS_DESCR_TYP_GROUPE'  pour un descripteur de type "groupe"
 *----------------------------------------------------------------------------*/

ecs_descr_typ_t
ecs_descr__ret_typ(const ecs_descr_t  *descr)
{
  assert(descr != NULL);

  return descr->typ;
}

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie l'identificateur du descripteur donné en argument
 *----------------------------------------------------------------------------*/

int
ecs_descr__ret_ide(const ecs_descr_t  *descr)
{
  assert(descr != NULL);

  return descr->ide;
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

