/*============================================================================
 *  Définitions des fonctions
 *   associées à la structure `ecs_table_t' décrivant une table
 *   et propres aux tables auxiliaires de type "attribut"
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
#include <string.h> /* strlen() */


/*----------------------------------------------------------------------------
 *  Fichiers `include' système
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/

#include "ecs_def.h"
#include "ecs_mem.h"
#include "ecs_tab.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles des paquetages visibles
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_table.h"
#include "ecs_descr_chaine.h"
#include "ecs_descr.h"
#include "ecs_famille_chaine.h"


/*----------------------------------------------------------------------------
 *  Fichier  `include' du  paquetage courant associé au fichier courant
 *----------------------------------------------------------------------------*/

#include "ecs_table_att.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' prives   du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_table_priv.h"


/*============================================================================
 *                              Fonctions privées
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction qui renumérote à partir de `1' les références
 *   des numéros de descripteurs
 *   dans le tableau des valeurs d'un attribut donné
 *
 *  La fonction renvoie un tableau de renumérotation des descripteurs
 *   d'un attribut de référence donné
 *   Si le numéro de descripteur n'a pas été renumérote, on lui attribue
 *    la valeur `num_descr_defaut' dans le tableau de renumérotation
 *----------------------------------------------------------------------------*/

static ecs_tab_int_t
_table_att__renum_descr(ecs_table_t      *table_att,
                        const ecs_int_t   num_descr_defaut,
                        size_t            nbr_descr_old)
{
  size_t          cpt_descr;
  size_t          idescr_old;
  size_t          ival;
  size_t          nbr_val;
  ecs_int_t       num_descr_old;
  ecs_int_t     * renum_descr_old_new;

  ecs_tab_int_t   vect_transf_num_descr_new_old;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(table_att != NULL);

  ECS_MALLOC(renum_descr_old_new, nbr_descr_old, ecs_int_t);
  ECS_MALLOC(vect_transf_num_descr_new_old.val, nbr_descr_old, ecs_int_t);
  vect_transf_num_descr_new_old.nbr = nbr_descr_old;

  for (idescr_old = 0; idescr_old < nbr_descr_old; idescr_old++)
    renum_descr_old_new[idescr_old] = num_descr_defaut;

  nbr_val = ecs_table__ret_val_nbr(table_att);

  cpt_descr = 0;

  for (ival = 0; ival < nbr_val; ival++) {

    num_descr_old = table_att->val[ival] - 1;

    if (renum_descr_old_new[num_descr_old] == num_descr_defaut) {

      renum_descr_old_new[num_descr_old] = cpt_descr;
      table_att->val[ival] = cpt_descr + 1;
      vect_transf_num_descr_new_old.val[cpt_descr++] = num_descr_old;

    }
    else
      table_att->val[ival] = renum_descr_old_new[num_descr_old] + 1;

  }

  ECS_FREE(renum_descr_old_new);

  ECS_REALLOC(vect_transf_num_descr_new_old.val, cpt_descr, ecs_int_t);
  vect_transf_num_descr_new_old.nbr = cpt_descr;

  return vect_transf_num_descr_new_old;
}

/*----------------------------------------------------------------------------
 *  Fonction réalisant la mise à jour des valeurs d'un attribut des élément
 *   lorque ceux-ci sont renumérotés. Chaque élément d'origine a au
 *   maximum un élément correspondant, mais un nouvel élément peut résulter
 *   de la fusion de plusieurs éléments d'origine, et en conserver tous
 *   les attributs.
 *----------------------------------------------------------------------------*/

static void
_table_att__herite(ecs_table_t          *table_att_elt,
                   const ecs_tab_int_t  *tab_old_new)
{
  int           ival_elt;
  ecs_int_t     ival_new;
  size_t        ielt;

  size_t        nbr_val_new;
  size_t        nbr_new;

  ecs_size_t    ival_att_new;
  ecs_int_t     ielt_new;
  ecs_int_t     num_elt_new;

  ecs_size_t    pos_val_att;

  ecs_int_t     val_elt;

  int         *cpt_att_new = NULL;
  ecs_size_t  *pos_att_new = NULL;
  ecs_int_t   *val_att_new = NULL;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(table_att_elt != NULL);

  /* Count values */

  ECS_MALLOC(pos_att_new, tab_old_new->nbr + 1, ecs_size_t);

  nbr_new = 0;

  for (ielt = 0; ielt < tab_old_new->nbr; ielt++)
    pos_att_new[ielt] = 0;

  for (ielt = 0; ielt < tab_old_new->nbr; ielt++) {
    num_elt_new = tab_old_new->val[ielt];
    if (num_elt_new != 0) {
      num_elt_new = ECS_ABS(num_elt_new); /* 1 to n */
      if (num_elt_new > (ecs_int_t)nbr_new)
        nbr_new = num_elt_new;
      pos_att_new[num_elt_new] +=   table_att_elt->pos[ielt+1]
                                  - table_att_elt->pos[ielt];
    }
  }

  /* Transform to position array (index) */

  pos_att_new[0] = 1;
  for (ielt = 0; ielt < nbr_new; ielt++)
    pos_att_new[ielt+1] += pos_att_new[ielt];

  nbr_val_new = pos_att_new[nbr_new] - 1;

  /* Now copy values */

  if (nbr_val_new != 0) {

    ECS_MALLOC(val_att_new, nbr_val_new, ecs_int_t);

    ECS_MALLOC(cpt_att_new, nbr_new, int);
    for (ielt = 0; ielt < nbr_new; ielt++)
      cpt_att_new[ielt] = 0;

    for (ielt = 0; ielt < tab_old_new->nbr; ielt++) {

      num_elt_new = tab_old_new->val[ielt];

      if (num_elt_new != 0) {

        ielt_new = ECS_ABS(num_elt_new) - 1;
        ival_att_new = pos_att_new[ielt_new] - 1;

        for (pos_val_att = table_att_elt->pos[ielt];
             pos_val_att < table_att_elt->pos[ielt+1];
             pos_val_att++) {

          /* Value to add */

          val_elt = table_att_elt->val[pos_val_att - 1];

          /* Add it only if it is not already present */

          for (ival_new = 0;
               (   ival_new < cpt_att_new[ielt_new]
                && val_att_new[ival_att_new + ival_new] != val_elt);
               ival_new++);

          if (ival_new == cpt_att_new[ielt_new]) {
            val_att_new[ival_att_new + ival_new] = val_elt;
            cpt_att_new[ielt_new] += 1;
          }
        }
      } /* End of addition if num_elt_new != 0) */

    } /* End of loop on old elements */

    /* Compact definitions */

    ival_new = 0;
    for (ielt = 0; ielt < nbr_new; ielt++) {

      ival_att_new = pos_att_new[ielt] - 1;

      for (ival_elt = 0; ival_elt < cpt_att_new[ielt]; ival_elt++)
        val_att_new[ival_new++] = val_att_new[ival_att_new + ival_elt];
    }

    for (ielt = 0; ielt < nbr_new; ielt++)
      pos_att_new[ielt + 1] = pos_att_new[ielt] + cpt_att_new[ielt];

    nbr_val_new = pos_att_new[nbr_new] - 1;

    ECS_FREE(cpt_att_new);
  }

  /* Update definitions */

  ECS_FREE(table_att_elt->pos);
  ECS_FREE(table_att_elt->val);

  ECS_REALLOC(pos_att_new, nbr_new + 1, ecs_size_t);
  ECS_REALLOC(val_att_new, nbr_val_new, ecs_int_t);

  table_att_elt->nbr = nbr_new;
  table_att_elt->pos = pos_att_new;
  table_att_elt->val = val_att_new;
}

/*----------------------------------------------------------------------------
 *  Fonction qui remplace les références à des éléments
 *  en des références à d'autres éléments liés aux premiers
 *  par un vecteur indexé contenant des valeurs positives mais pouvant
 *  contenir la valeur 0
 *----------------------------------------------------------------------------*/

static ecs_table_t *
_table_att__famille_en_groupe(size_t               n_elts,
                              const int           *elt_fam,
                              const ecs_table_t   *table_def)
{
  ecs_table_t  *table_rep_new;

  size_t        cpt_val_rep;
  size_t        ielt_rep;
  size_t        ipos_def;
  size_t        nbr_elt_rep;
  size_t        nbr_pos_def;
  size_t        nbr_val_rep_new;
  ecs_int_t     num_val_def;
  size_t        pos_deb_def;
  size_t        pos_fin_def;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  nbr_elt_rep = n_elts;

  /* Passage pour dimensionnement */

  cpt_val_rep = 0;
  nbr_val_rep_new = 0;

  for (ielt_rep = 0; ielt_rep < nbr_elt_rep; ielt_rep++) {

    num_val_def = elt_fam[ielt_rep];

    if (num_val_def > 0) {

      pos_deb_def = table_def->pos[num_val_def - 1] - 1;
      pos_fin_def = table_def->pos[num_val_def    ] - 1;

      nbr_pos_def = pos_fin_def - pos_deb_def;
    }

    else if (num_val_def == 0)
      nbr_pos_def = 0;

    else
      assert(num_val_def >= 0);

    nbr_val_rep_new = ECS_MAX(nbr_val_rep_new, cpt_val_rep + nbr_pos_def);

    cpt_val_rep += nbr_pos_def;

  }

  /* Traitement effectif */

  table_rep_new = ecs_table__alloue(nbr_elt_rep,
                                    nbr_val_rep_new);

  table_rep_new->pos[0] = 1;

  cpt_val_rep = 0;

  for (ielt_rep = 0; ielt_rep < nbr_elt_rep; ielt_rep++) {

    num_val_def = elt_fam[ielt_rep];

    if (num_val_def > 0) {

      pos_deb_def = table_def->pos[num_val_def - 1] - 1;
      pos_fin_def = table_def->pos[num_val_def    ] - 1;

      nbr_pos_def = pos_fin_def - pos_deb_def;

      for (ipos_def = 0; ipos_def < nbr_pos_def; ipos_def++)
        table_rep_new->val[cpt_val_rep + ipos_def]
          = table_def->val[pos_deb_def + ipos_def];

    }
    else if (num_val_def == 0) {
      nbr_pos_def = 0;
    }
    else {
      assert(num_val_def >= 0);
    }

    cpt_val_rep += nbr_pos_def;

    table_rep_new->pos[ielt_rep + 1] = cpt_val_rep + 1;
  }

  return table_rep_new;
}

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie un vecteur indexé donnant pour chaque élément
 *   le numéro de famille à laquelle il appartient
 *   à partir -    du vecteur indexé donnant pour chaque élément
 *               les numéros d'attribut le caractérisant
 *            - et du tableau donnant pour chaque élément :
 *                 - le 1er élément ayant la même liste d'attributs
 *                 - ou `-1' s'il n'a pas d'attribut
 *
 *  La fonction détermine aussi :
 *   - les définitions des familles en fonction des numéros de descripteur
 *   - les nombres de descripteurs composant chaque famille
 *   - le  nombre  de familles
 *----------------------------------------------------------------------------*/

static int *
_table__attribue_fam(const ecs_table_t      *table_att_unifie,
                     const ecs_tab_int_t     tab_idem,
                     ecs_int_t            ***def_fam_descr,
                     ecs_int_t             **nbr_descr_fam,
                     const ecs_int_t         num_fam_deb,
                     int                    *nbr_fam)
{
  ecs_int_t     cpt_fam;
  ecs_int_t     idescr;
  size_t        ielt;
  ecs_int_t     ifam;
  ecs_int_t     nbr_att_elt;
  size_t        nbr_elt;
  size_t        pos_elt;

  int *elt_fam;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  nbr_elt = table_att_unifie->nbr;

  ECS_MALLOC(elt_fam, nbr_elt, int);

  /* Détermination du nombre de familles */
  /*-------------------------------------*/

  *nbr_fam = 0;

  for (ielt = 0; ielt < nbr_elt; ielt++) {
    if (tab_idem.val[ielt] == (ecs_int_t)ielt)
      (*nbr_fam)++;
  }

  /* Détermination du nombre de descripteurs par famille */
  /*-----------------------------------------------------*/

  ECS_MALLOC(*nbr_descr_fam, *nbr_fam, ecs_int_t);

  cpt_fam = 0;

  for (ielt = 0; ielt < nbr_elt; ielt++) {

    if (tab_idem.val[ielt] == (ecs_int_t)ielt) {

      pos_elt = table_att_unifie->pos[ielt] - 1;
      (*nbr_descr_fam)[cpt_fam]
        = table_att_unifie->pos[ielt + 1] - 1 - pos_elt;

      cpt_fam++;
    }
  }

  /* Détermination des définitions des familles en fonction des descripteurs */
  /*-------------------------------------------------------------------------*/

  ECS_MALLOC(*def_fam_descr, *nbr_fam, ecs_int_t *);

  for (ifam = 0; ifam < *nbr_fam; ifam++)
    ECS_MALLOC((*def_fam_descr)[ifam], (*nbr_descr_fam)[ifam], ecs_int_t);

  cpt_fam = 0;

  for (ielt = 0; ielt < nbr_elt; ielt++) {

    if (tab_idem.val[ielt] == (ecs_int_t)ielt) {

      /* Définition de la famille en fonction des numéros d'attribut */

      pos_elt     = table_att_unifie->pos[ielt    ] - 1;
      nbr_att_elt = table_att_unifie->pos[ielt + 1] - 1 - pos_elt;

      for (idescr = 0; idescr < nbr_att_elt; idescr++)
        (*def_fam_descr)[cpt_fam][idescr]
          = table_att_unifie->val[pos_elt + idescr];

      cpt_fam++;
    }

    if (tab_idem.val[ielt] != -1) {

      elt_fam[ielt] = cpt_fam + num_fam_deb - 1;

    }
    else {

      /* L'élément n'a pas de famille */

      elt_fam[ielt] = 0;
    }
  }

  return elt_fam;
}

/*----------------------------------------------------------------------------
 *  Fonction qui compare les ensembles de valeurs comprises entre 2 positions
 *   et renvoie un tableau donnant pour chaque ensemble
 *      le numéro du 1er ensemble de chaque groupe d'ensembles identiques
 *
 *  Si entre 2 positions il n'y a pas de valeurs, le numéro est `-1'
 *
 *  Cette fonction prend tout son intérêt lorsque
 *   - les valeurs ont été triées dans chaque ensemble de valeurs
 *   - les ensembles de valeurs ont été tries dans l'ordre lexicographique
 *----------------------------------------------------------------------------*/

static ecs_tab_int_t
_table_att__compare_val_pos(const ecs_table_t  *this_table)
{
  bool          bool_diff;

  size_t        ipos;
  size_t        ival;
  size_t        max_val_pos;
  size_t        nbr_val_pos;
  size_t        nbr_val_pos_ref;
  size_t        pos_ref;
  size_t        pos_val;

  ecs_int_t   * val_pos_ref;

  ecs_tab_int_t   tab_idem;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /*=================*/
  /* Initialisations */
  /*=================*/

  tab_idem.nbr = this_table->nbr;
  ECS_MALLOC(tab_idem.val, tab_idem.nbr, ecs_int_t);

  max_val_pos = this_table->pos[1] - this_table->pos[0];
  for (ipos = 1; ipos < tab_idem.nbr; ipos++)
    max_val_pos = ECS_MAX(max_val_pos,
                          this_table->pos[ipos + 1] -
                          this_table->pos[ipos    ]   );

  ECS_MALLOC(val_pos_ref, max_val_pos, ecs_int_t);

  /*=================================================================*/
  /* Repérage des ensembles de valeurs (comprises entre 2 positions) */
  /*  identiques par comparaison de deux ensembles successifs        */
  /*=================================================================*/

  /*----------------------------------------------*/
  /* Initialisation pour le tout premier ensemble */
  /*----------------------------------------------*/

  /* Ensemble servant de référence au 1er tour */

  pos_val         = this_table->pos[0] - 1;
  nbr_val_pos_ref = this_table->pos[1] - 1 - pos_val;

  for (ival = 0; ival < nbr_val_pos_ref; ival++)
    val_pos_ref[ival] = this_table->val[pos_val + ival];

  pos_ref = 0;

  if (nbr_val_pos_ref == 0)
    tab_idem.val[0] = -1;
  else
    tab_idem.val[0] =  0;

  /*---------------------------------*/
  /* Boucle sur les autres ensembles */
  /*---------------------------------*/

  for (ipos = 1; ipos < tab_idem.nbr; ipos++) {

    pos_val     = this_table->pos[ipos    ] - 1;
    nbr_val_pos = this_table->pos[ipos + 1] - 1 - pos_val;

    if (nbr_val_pos == nbr_val_pos_ref) {

      ival = 0;
      while (ival < nbr_val_pos                                &&
             ECS_ABS(this_table->val[pos_val + ival]) ==
             ECS_ABS(val_pos_ref[ival])                           )
        ival++;

      if (ival == nbr_val_pos) {

        /* Les 2 ensembles de valeurs consécutifs sont identiques */

        tab_idem.val[ipos] = tab_idem.val[pos_ref];

        bool_diff = false;

      }
      else {

        /* Les 2 ensembles de valeurs consécutifs sont différents */

        bool_diff = true;

      }

    }
    else{

      /* Les 2 ensembles de valeurs consécutifs sont différents */

      bool_diff = true;

    }

    if (bool_diff == true) { /* Si les deux ensembles sont différents */

      /* L'ensemble qui vient d'être comparé       */
      /*  devient l'ensemble de référence          */
      /*  pour la comparaison du prochain ensemble */

      nbr_val_pos_ref = nbr_val_pos;

      for (ival = 0; ival < nbr_val_pos; ival++)
        val_pos_ref[ival] = this_table->val[pos_val + ival];

      pos_ref = ipos;

      if (nbr_val_pos == 0)
        tab_idem.val[ipos] = -1;
      else
        tab_idem.val[ipos] = ipos;

    } /* Fin : si les 2 ensembles sont différents */

  } /* Fin : boucle sur les autres ensembles */

  ECS_FREE(val_pos_ref);

  return tab_idem;
}

/*----------------------------------------------------------------------------
 *    Fonction de descente d'un arbre binaire pour le tri lexicographique
 *  d'une table de listes d'entiers (voir ecs_int_table_ord).
 *----------------------------------------------------------------------------*/

static void
_table_att__desc_arbre(const ecs_size_t  *pos_tab,
                       const ecs_int_t   *val_tab,
                       size_t             ltree,
                       size_t             ntree,
                       ecs_int_t         *elem_ord)
{
  size_t    ktree;
  size_t    p1, p2, l1, l2, i, i_save;
  size_t    e1, e2;
  int     isup;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  i_save = *(elem_ord+ltree);

  while (ltree <= (ntree/2)) {

    ktree = (2*ltree)+1;

    if (ktree < ntree - 1) {

      /* isup vrai si elem_ord[ktree+1] > elem_ord[ktree] */
      p1 = *(elem_ord+ktree+1);
      p2 = *(elem_ord+ktree);
      l1 = *(pos_tab+p1+1) - *(pos_tab+p1);
      l2 = *(pos_tab+p2+1) - *(pos_tab+p2);
      isup = (l1 > l2) ? 1 : 0;
      for (i = 0, e1 = *(pos_tab+p1) - 1, e2 = *(pos_tab+p2) - 1;
           (i < l1) && (i < l2)
             && (*(val_tab+e1) == *(val_tab+e2));
           i++, e1++, e2++);
      if ((i < l1) && (i < l2))
        isup = (*(val_tab+e1) > *(val_tab+e2)) ? 1 : 0;

      /* si isup vrai, on incremente ktree */
      if (isup) ktree++;
    }

    if (ktree >= ntree) break;

    /* isup faux si elem_ord[ltree] (initial) < elem_ord[ktree] */
    p1 = i_save;
    p2 = *(elem_ord+ktree);
    l1 = *(pos_tab+p1+1) - *(pos_tab+p1);
    l2 = *(pos_tab+p2+1) - *(pos_tab+p2);
    isup = (l1 < l2) ? 0 : 1;
    for (i = 0, e1 = *(pos_tab+p1) - 1, e2 = *(pos_tab+p2) - 1;
         (i < l1) && (i < l2) &&
           (*(val_tab+e1) == *(val_tab+e2));
         i++, e1++, e2++);
    if ((i < l1) && (i < l2))
      isup = (*(val_tab+e1) < *(val_tab+e2)) ? 0 : 1;

    /* si isup vrai */
    if (isup) break;

    *(elem_ord+ltree) = *(elem_ord+ktree);
    ltree = ktree;
  }

  *(elem_ord+ltree) = i_save;
}

/*----------------------------------------------------------------------------
 *    Fonction de tri lexicographique d'une table de listes d'entiers.
 *  La liste n'est pas modifiée directement, mais on construit une table de
 *  renumérotation, afin de pouvoir appliquer cette renumérotation à d'autres
 *  tableaux dépendant de la table.
 *
 *    Le tri utilisé est de type "heapsort", de complexité O(nlog(n)). Les
 *  éléments sont rangés en ordre croissant.
 *
 *    Pour accéder au k-ième élément ordonné de la table,
 *  on utilisera la position :
 *              elem_pos[elem_ord[k - 1] - 1] du tableau elem_table[].
 *
 * ----------------------------------------------------------------------------
 *
 *  Exemple : (pour des tables de valeurs entières)
 *  =======
 *
 *  Supposons 7 éléments définis par des sommets de numéros 1,2,3,4,5,6,9
 *   comme suit :
 *
 *                        .---.---..---..---.---.---..---..---.---..---..---.
 *  this_vec->val         | 5 | 3 || 4 || 5 | 2 | 6 || 4 || 5 | 9 || 1 || 2 |
 *                        `---'---'`---'`---'---'---'`---'`---'---'`---'`---'
 *                          0   1    2    3   4   5    6    7   8    9    10
 *                        `---.---'`-.-'`-----.-----'`-.-'`---.---'`-.-'`-.-'
 *  Numéro de l'élément :     1      2        3        4      5      6    7
 *
 *
 *                        .---.---.---.---.---.---.---.---.
 *  this_vec->pos         | 1 | 3 | 4 | 7 | 8 | 10| 11| 12|
 *                        `---'---'---'---'---'---'---'---'
 *                          0   1   2   3   4   5   6   7
 *
 *
 *
 *  On veut trier les 6 éléments 1,2,3,4,5,7
 *  ----------------------------------------
 *
 *
 *                        .---.---.---.---.---.---.
 *  vect_renum            | 0 | 1 | 2 | 3 | 4 | 6 |
 *                        `---'---'---'---'---'---'
 *                          0   1   2   3   4   5
 *
 *
 *  La liste trieé étant :
 *  --------------------
 *                        .---..---..---..---..---.---.---..---.---..---.---.
 *                        | 1 || 2 || 4 || 4 || 5 | 2 | 6 || 5 | 3 || 5 | 9 |
 *                        `---'`---'`---'`---'`---'---'---'`---'---'`---'---'
 *                        `-.-'`-.-'`-.-'`-.-'`-----.-----'`---.---'`---.---'
 *       Ancien  numéro :   6    7    2    4        3          1        5
 *       Nouveau numéro :   1    2    3    4        5          6        7
 *
 *
 *
 *  on obtient pour résultat :
 *  ------------------------
 *
 *
 *                        .---.---.---.---.---.---.
 *  vect_renum            | 6 | 1 | 3 | 2 | 0 | 4 |
 *                        `---'---'---'---'---'---'
 *                          0   1   2   3   4   5
 *
 *----------------------------------------------------------------------------*/

static void
_table_att__trie(const ecs_table_t  *this_table,
                 ecs_tab_int_t      *vect_renum)
{
  ecs_int_t i, i_save;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* création de l'arbre binaire vec_renum->val[vect_renum->nbr] */

  for (i = (vect_renum->nbr / 2) - 1; i >= 0; i--) {

    _table_att__desc_arbre(this_table->pos,
                           this_table->val,
                           i,
                           vect_renum->nbr,
                           vect_renum->val);
  }

  /* tri de l'arbre binaire */

  for (i = vect_renum->nbr - 1; i > 0; i--) {

    i_save                 = *(vect_renum->val    );
    *(vect_renum->val    ) = *(vect_renum->val + i);
    *(vect_renum->val + i) = i_save;

    _table_att__desc_arbre(this_table->pos,
                           this_table->val,
                           0,
                           i,
                           vect_renum->val);
  }
}

/*----------------------------------------------------------------------------
 *  Fonction qui trie un vecteur d'entiers en renvoyant le vecteur trié
 *
 *  --------------------------------------------------------------------------
 *
 *  Suite de l'exemple précédent :
 *  ============================
 *
 *  Cette fonction renvoie, en plus, le vecteur trié :
 *  ------------------------------------------------
 *
 *                        .---..---..---..---.---.---..---.---..---.---.
 *  vec_trie->val         | 2 || 4 || 4 || 5 | 2 | 6 || 5 | 3 || 5 | 9 |
 *                        `---'`---'`---'`---'---'---'`---'---'`---'---'
 *                          0    1    2    3   4   5    6   7    8   9
 *                        `-,-'`-,-'`-,-'`-----,-----'`---,---'`---,---'
 *       Ancien  numéro :   7    2    4        3          1        5
 *       Nouveau numéro :   1    2    3        4          5        6
 *
 *
 *                        .---.---.---.---.---.---.---.
 *  vec_trie->pos         | 1 | 2 | 3 | 4 | 7 | 9 | 11|
 *                        `---'---'---'---'---'---'---'
 *                          0   1   2   3   4   5   6
 *
 *----------------------------------------------------------------------------*/

static void
_table_att__trie_et_renvoie(const ecs_table_t  *this_table,
                            ecs_table_t        *table_trie,
                            ecs_tab_int_t      *vect_renum_pos)
{
  size_t      pos_nbr_val;
  size_t      ipos;
  size_t      ival;
  size_t      ival_deb;
  size_t      cpt_val;
  ecs_int_t * vect_renum_pos_val;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  _table_att__trie(this_table, vect_renum_pos);

  cpt_val = 0;

  table_trie->nbr = vect_renum_pos->nbr;

  table_trie->pos[0] = 1;

  for (ipos = 0; ipos < vect_renum_pos->nbr; ipos++) {

    pos_nbr_val
      = this_table->pos[vect_renum_pos->val[ipos] + 1]
      - this_table->pos[vect_renum_pos->val[ipos]];

    table_trie->pos[ipos + 1] = table_trie->pos[ipos] + pos_nbr_val;

    ival_deb = this_table->pos[vect_renum_pos->val[ipos]] - 1;

    for (ival = 0; ival < pos_nbr_val; ival++)
      table_trie->val[cpt_val++] = this_table->val[ival_deb+ival];
  }

  /* `vect_renum_pos' prend pour indice les indices nouveaux,       */
  /*  et ses valeurs contiennent les indices anciens correspondants */
  /* On inverse le contenu de `vect_renum_pos' :                    */
  /*  à chaque indice ancien, `vect_renum_pos' donne la valeur      */
  /*  du nouvel indice                                              */

  ECS_MALLOC(vect_renum_pos_val, vect_renum_pos->nbr, ecs_int_t);

  for (ival = 0; ival < vect_renum_pos->nbr; ival++)
    vect_renum_pos_val[ival] = vect_renum_pos->val[ival];

  for (ival = 0; ival < vect_renum_pos->nbr; ival++)
    vect_renum_pos->val[vect_renum_pos_val[ival]] = ival;

  ECS_FREE(vect_renum_pos_val);
}

/*----------------------------------------------------------------------------
 *    Fonction de descente d'un arbre binaire pour le tri sur place
 *  d'un tableau d'entiers.
 *----------------------------------------------------------------------------*/

static void
_table_att__desc_arbre_val(ecs_int_t  ltree,
                           ecs_int_t  ntree,
                           ecs_int_t  elem_ord[])
{
  ecs_int_t ktree;
  ecs_int_t i_save;

  i_save = elem_ord[ltree];

  while (ltree <= (ntree / 2)) {

    ktree = (2 * ltree) + 1;

    if (ktree < ntree - 1)
      if (elem_ord[ktree + 1] > elem_ord[ktree]) ktree++;

    if (ktree >= ntree) break;

    if (i_save >= elem_ord[ktree]) break;

    elem_ord[ltree] = elem_ord[ktree];
    ltree = ktree;

  }

  elem_ord[ltree] = i_save;
}

/*----------------------------------------------------------------------------
 *  Fonction qui réalise un tri sur les valeurs d'un vecteur indexé
 *   entre 2 positions
 *----------------------------------------------------------------------------*/

static ecs_table_t  *
_table_att__trie_val_pos(const ecs_table_t  *this_table)
{
  ecs_table_t *table_ord;

  ecs_int_t    nbr_elt;
  ecs_int_t    nbr_def;
  ecs_int_t    ielt;
  ecs_int_t    ipos;
  ecs_int_t    ipos_deb;
  ecs_int_t    ipos_fin;
  ecs_int_t    nbr_def_loc;
  ecs_int_t    i;
  ecs_int_t    i_save;

  ecs_int_t   *tab_val_loc;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /*=================*/
  /* Initialisations */
  /*=================*/

  nbr_elt = this_table->nbr;
  nbr_def = this_table->pos[nbr_elt] - 1;

  table_ord = ecs_table__alloue(nbr_elt,
                                nbr_def);

  table_ord->pos[0] = 1;

  /*=========================*/
  /* Boucle sur les éléments */
  /*=========================*/

  for (ielt = 0; ielt < nbr_elt; ielt++) {

    /* Les positions restent inchangées */

    table_ord->pos[ielt + 1] = this_table->pos[ielt + 1];

    ipos_deb = this_table->pos[ielt    ] - 1;
    ipos_fin = this_table->pos[ielt + 1] - 1;

    /* On commence par recopier les valeurs */

    for (ipos = ipos_deb; ipos < ipos_fin; ipos++)
      table_ord->val[ipos] = this_table->val[ipos];

    /* Tri des valeurs de la copie */
    /*-----------------------------*/

    nbr_def_loc = ipos_fin -ipos_deb;

    tab_val_loc = table_ord->val + ipos_deb;

    /* Création de l'arbre binaire table_renum->val[vect_renum->nbr] */

    for (i = (nbr_def_loc / 2) - 1; i >= 0; i--)
      _table_att__desc_arbre_val(i, nbr_def_loc, tab_val_loc);

    /* Tri de l'arbre binaire */

    for (i = nbr_def_loc - 1; i > 0; i--) {

      i_save         = tab_val_loc[0];
      tab_val_loc[0] = tab_val_loc[i];
      tab_val_loc[i] = i_save;

      _table_att__desc_arbre_val(0, i, tab_val_loc);
    }

  } /* Fin : boucle sur les éléments */

  /* Renvoi de la structure des définitions ordonnées */
  /*--------------------------------------------------*/

  return table_ord;
}

/*----------------------------------------------------------------------------
 *  Fonction qui réalise l'ordination des valeurs d'un vecteur indexé
 *     en triant les valeurs pour chaque ensemble de valeurs entre 2 positions
 *  et en triant chaque ensemble de valeurs dans l'ordre lexicographique
 *----------------------------------------------------------------------------*/

static ecs_tab_int_t
_table_att__trie_val(ecs_table_t  **this_table)
{
  ecs_table_t  *table_def_ord;
  ecs_table_t  *table_elt_ord;
  ecs_tab_int_t   vect_transf;
  size_t   nbr_elt;
  size_t   nbr_val;
  size_t   ielt;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  nbr_elt = (*this_table)->nbr;

  /* Tri sur les définitions des éléments */

  table_def_ord = _table_att__trie_val_pos(*this_table);

  /* Tri des éléments en fonction de leurs définitions ordonnées */

  nbr_val = ecs_table__ret_val_nbr(table_def_ord);

  vect_transf.nbr = nbr_elt;
  ECS_MALLOC(vect_transf.val, nbr_elt, ecs_int_t);
  for (ielt = 0; ielt < nbr_elt; ielt++)
    vect_transf.val[ielt] = ielt;

  table_elt_ord = ecs_table__alloue(nbr_elt,  nbr_val);

  table_elt_ord->descr = (*this_table)->descr;
  (*this_table)->descr = NULL;

  ecs_table__detruit(this_table);

  _table_att__trie_et_renvoie(table_def_ord,
                              table_elt_ord,
                              &vect_transf);

  ecs_table__detruit(&table_def_ord);

  /* Affectation de la nouvelle liste compactée */

  *this_table = table_elt_ord;

  /* Renvoi du vecteur de transformation */

  return vect_transf;
}

/*----------------------------------------------------------------------------
 *  Fonction qui trie les descripteurs d'une table de type "attribut"
 *   et met à jour le numérotation correspondante des valeurs.
 *----------------------------------------------------------------------------*/

static void
_table_att__trie_descr(ecs_table_t  *table_att)
{
  ecs_int_t      nbr_val, ival;
  ecs_tab_int_t  tab_renum_descr;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(table_att != NULL);

  tab_renum_descr = ecs_descr_chaine__trie(table_att->descr);

  nbr_val = ecs_table__ret_val_nbr(table_att);

  for (ival = 0; ival < nbr_val; ival++) {
    ecs_int_t idescr = table_att->val[ival] -1;
    table_att->val[ival] = tab_renum_descr.val[idescr];
  }

  tab_renum_descr.nbr = 0;
  ECS_FREE(tab_renum_descr.val);
}

/*----------------------------------------------------------------------------
 *  Fonction qui retourne le tableau des valeurs d'un tableau donné,
 *   dimensionné au nombre de valeurs distinctes de ce tableau.
 *----------------------------------------------------------------------------*/

static ecs_tab_int_t
_table_att__ret_reference(size_t   n_elts,
                          int     *elt_fam)
{
  size_t        cpt_ref;
  size_t        nbr_ref;
  ecs_int_t     ind_ref;
  size_t        iref;
  size_t        ival;
  ecs_int_t    *liste_ref;
  ecs_tab_int_t   tab_ref;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Dimensionnement a priori de le liste des références */

  nbr_ref = 10;
  ECS_MALLOC(liste_ref, nbr_ref, ecs_int_t);
  for (iref = 0; iref < nbr_ref; iref++)
    liste_ref[iref] = false;

  for (ival = 0; ival < n_elts; ival++) {

    ind_ref = elt_fam[ival] - 1;

    /* Re-dimensionnement si nécessaire de la liste des références */

    if (ind_ref >= (ecs_int_t)nbr_ref) {
      ECS_REALLOC(liste_ref, ind_ref * 2, ecs_int_t);
      for (iref = nbr_ref; iref < (size_t)ind_ref*2; iref++)
        liste_ref[iref] = false;
      nbr_ref = ind_ref * 2;
    }

    if (ind_ref > -1)
      liste_ref[ind_ref] = true;

  }

  ECS_MALLOC(tab_ref.val, nbr_ref, ecs_int_t);

  cpt_ref = 0;

  for (iref = 0; iref < nbr_ref; iref++) {
    if (liste_ref[iref] == true)
      tab_ref.val[cpt_ref++] = iref + 1;
  }

  ECS_FREE(liste_ref);

  if (cpt_ref != 0)
    ECS_REALLOC(tab_ref.val, cpt_ref, ecs_int_t);
  else
    ECS_FREE(tab_ref.val);
  tab_ref.nbr = cpt_ref;

  return tab_ref;
}

/*----------------------------------------------------------------------------
 *  Fonction qui créé une structure `ecs_table_t'
 *   à partir d'un tableau bi-dimensionnel `tab_elt' contenant
 *   pour chaque élément un ensemble de valeurs
 *  Si un élément `ielt' n'a pas de valeur associée, `tab_elt[ielt].nbr = 0'
 *----------------------------------------------------------------------------*/

static ecs_table_t *
_table_att__transforme_bi_tab(const ecs_tab_int_t  *tab_elt,
                              size_t                nbr_elt,
                              size_t                nbr_val)
{
  size_t        cpt_val;
  size_t        ielt;
  size_t        ival;

  ecs_table_t * this_table;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(tab_elt != NULL);

  /* Construction des tableaux de positions et de valeurs */
  /*------------------------------------------------------*/

  this_table = ecs_table__alloue(nbr_elt, nbr_val);

  this_table->pos[0] = 1;
  cpt_val = 0;

  for (ielt = 0; ielt < nbr_elt; ielt++) {

    if (tab_elt[ielt].nbr != 0) {

      this_table->pos[ielt + 1]
        = this_table->pos[ielt] + tab_elt[ielt].nbr;

      for (ival = 0; ival < tab_elt[ielt].nbr; ival++)
        this_table->val[cpt_val++] = tab_elt[ielt].val[ival];

    }
    else
      this_table->pos[ielt + 1]  = this_table->pos[ielt];
  }

  return this_table;
}

/*----------------------------------------------------------------------------
 *  Fonction réalisant la transformation d'un tableau de valeurs de familles
 *   en appliquant directement le vecteur de transformation donné
 *
 *  Le nombre de valeurs transformées doit être égal
 *  au nombre de valeurs avant transformation
 *----------------------------------------------------------------------------*/

static void
_ecs_table__transforme_fam(size_t                n_elts,
                           int                  *elt_fam,
                           const ecs_tab_int_t   vect_transf)
{
  int         *fam_ref;
  size_t       ielt;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (elt_fam == NULL)
    return;

  assert(vect_transf.nbr == n_elts);

  ECS_MALLOC(fam_ref, n_elts, int);

  memcpy(fam_ref, elt_fam, n_elts*sizeof(int));

  for (ielt = 0; ielt < n_elts; ielt++) {

    assert(vect_transf.val[ielt] > -1);
    assert(vect_transf.val[ielt] < (ecs_int_t)n_elts);

    elt_fam[ielt] = fam_ref[vect_transf.val[ielt]];

  }

  ECS_FREE(fam_ref);
}

/*============================================================================
 *                             Fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction qui assemble une table donnée dans une table réceptrice donnée
 *
 *  L'assemblage consiste à :
 *  - regrouper sous la même position les valeurs des 2 tables
 *    (cela suppose donc que les 2 tables ont le même nombre de positions)
 *  - assembler les membres des descripteurs des 2 tables
 *    Les descripteurs des 2 tables peuvent être à `NULL'
 *    et si le descripteur de la table réceptrice est à `NULL',
 *    le descripteur de la table assemblée est celui de la table à assembler
 *
 *  La table à assembler est détruite apres assemblage
 * ----------------------------------------------------------------------------
 *
 *  Exemple :
 *  =======
 *
 *  Soit la table à assembler :
 *
 *                         .---.---..---..---.---.---.
 *     assemb->val         | 5 | 3 || 4 || 5 | 2 | 6 |
 *                         `---'---'`---'`---'---'---'
 *                           0   1    2    3   4   5
 *
 *
 *                         .---.---.---.---.---.
 *     assemb->pos         | 1 | 3 | 4 | 4 | 7 |
 *                         `---'---'---'---'---'
 *                           0   1   2   3   4
 *
 *
 *  dans la table réceptrice :
 *
 *                         .---..---..---.---..---.
 *     recept->val         | 4 || 5 || 6 | 6 || 1 |
 *                         `---'`---'`---'---'`---'
 *                           0    1    2   3    4
 *
 *
 *                         .---.---.---.---.---.
 *     recept->pos         | 1 | 2 | 3 | 5 | 6 |
 *                         `---'---'---'---'---'
 *                           0   1   2   3   4
 *
 *
 *  La table réceptrice devient :
 *
 *                         .---.---.---..---.---..---.---..---.---.---.---.
 *     recept->val         | 4 | 5 | 3 || 5 | 4 || 6 | 6 || 1 | 5 | 2 | 6 |
 *                         `---'---'---'`---'---'`---'---'`---'---'---'---'
 *                           0   1   2    3   4    5   6    7   8   9   10
 *
 *
 *                         .---.---.---.---.---.
 *     recept->pos         | 1 | 4 | 6 | 8 | 12|
 *                         `---'---'---'---'---'
 *                           0   1   2   3   4
 *
 *----------------------------------------------------------------------------*/

void
ecs_table_att__assemble(ecs_table_t  *table_recept,
                        ecs_table_t  *table_assemb)
{
  size_t    ipos;
  size_t    ival;
  size_t    cpt_val;

  ecs_size_t * recept_pos_tab_new;
  ecs_int_t  * recept_val_tab_new;

  size_t      recept_val_nbr;
  size_t      assemb_val_nbr;

  ecs_descr_t * descr_tete_loc;

  int incr_num_descr = 0;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(table_recept != NULL);
  assert(table_assemb != NULL);
  assert(table_recept->nbr == table_assemb->nbr);

  /* Récupération du nombre de descripteurs du premier table si lieu */

  incr_num_descr = ecs_table__ret_descr_nbr(table_recept);

  /* Assemblage des valeurs */

  ecs_table__regle_en_pos(table_recept);
  ecs_table__regle_en_pos(table_assemb);

  /* Allocation des nouveaux tableaux de la structure receptrice */
  /*-------------------------------------------------------------*/

  recept_val_nbr = table_recept->pos[table_recept->nbr] - 1;
  assemb_val_nbr = table_assemb->pos[table_assemb->nbr] - 1;

  ECS_MALLOC(recept_pos_tab_new,
             table_recept->nbr + 1,
             ecs_size_t);

  ECS_MALLOC(recept_val_tab_new,
             recept_val_nbr + assemb_val_nbr,
             ecs_int_t);

  /* Assemblage des tables de positions et de valeurs */
  /*--------------------------------------------------*/

  cpt_val = 0;
  recept_pos_tab_new[0] = 1;

  for (ipos = 0; ipos < table_recept->nbr; ipos++) {

    recept_pos_tab_new[ipos + 1]
      =   recept_pos_tab_new[ipos]
        + table_recept->pos[ipos + 1] - table_recept->pos[ipos]
        + table_assemb->pos[ipos + 1] - table_assemb->pos[ipos];

    for (ival = (table_recept->pos[ipos    ] - 1);
         ival < (table_recept->pos[ipos + 1] - 1);
         ival++) {

      recept_val_tab_new[cpt_val++] = table_recept->val[ival];

    }

    for (ival = (table_assemb->pos[ipos    ] - 1);
         ival < (table_assemb->pos[ipos + 1] - 1);
         ival++)
      recept_val_tab_new[cpt_val++]
        = table_assemb->val[ival] + incr_num_descr;

  }

  ECS_FREE(table_recept->pos);
  ECS_FREE(table_recept->val);

  table_recept->pos = recept_pos_tab_new;
  table_recept->val = recept_val_tab_new;

  ecs_table__pos_en_regle(table_recept);

  /* Ajout des descripteurs de la table à assembler à la table réceptrice */

  descr_tete_loc = ecs_descr_chaine__copie(table_assemb->descr);

  ecs_descr_chaine__ajoute(&table_recept->descr,
                           descr_tete_loc);

  ecs_table__detruit(&table_assemb);
}

/*----------------------------------------------------------------------------
 *  Fonction réalisant la mise à jour des valeurs d'un attribut des élément
 *   lorque ceux-ci sont renumérotés. Chaque élément d'origine a au
 *   maximum un élément correspondant, mais un nouvel élément peut résulter
 *   de la fusion de plusieurs éléments d'origine, et en conserver tous
 *   les attributs.
 *----------------------------------------------------------------------------*/

void
ecs_table_att__herite(ecs_table_t    *table_att_elt,
                      ecs_tab_int_t  *tab_old_new)
{
  ecs_tab_int_t   vect_transf_num_descr;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (table_att_elt == NULL)
    return;

  ecs_table__regle_en_pos(table_att_elt);

  /* Construction, pour les éléments,                    */
  /*  du vecteur `ecs_table_t' associé au table attribut */
  /*-----------------------------------------------------*/

  _table_att__herite(table_att_elt, tab_old_new);

  /* Allocation et initialisation pour les sous-éléments,   */
  /*  des vecteurs `ecs_table_t' associé au table attribut  */
  /*--------------------------------------------------------*/

  if (table_att_elt->descr != NULL) {

    const ecs_int_t nbr_descr_old
      = ecs_descr_chaine__ret_nbr(table_att_elt->descr);

    ecs_descr_t *descr_old = table_att_elt->descr;

    /* Renumérotation des numéros de descripteur */
    /*-------------------------------------------*/

    vect_transf_num_descr = _table_att__renum_descr(table_att_elt,
                                                    ECS_DESCR_NUM_NUL,
                                                    nbr_descr_old);

    table_att_elt->descr
      = ecs_descr_chaine__renumerote(table_att_elt->descr,
                                     vect_transf_num_descr );

    ecs_descr_chaine__detruit(&descr_old);

    ECS_FREE(vect_transf_num_descr.val);
  }

  ecs_table__pos_en_regle(table_att_elt);
}

/*----------------------------------------------------------------------------
 *  Fonction réalisant la transformation d'une table
 *   en fusionnant les propriétés de ses éléments
 *   qui sont identiquement transformes par le vecteur de transformation donné
 *----------------------------------------------------------------------------*/

void
ecs_table_att__fusionne(ecs_table_t          *this_table_att,
                        size_t                nbr_elt_new,
                        const ecs_tab_int_t   vect_transf)
{
  ecs_size_t  * pos_tab_transf;
  ecs_int_t  ** val_tab_transf;
  ecs_int_t   * val_tab_unidim;

  size_t       cpt_val_transf;
  size_t       nbr_elt_ref;
  size_t       nbr_val_transf;
  ecs_int_t    num_elt_transf;
  size_t       pos_inf;
  size_t       pos_sup;
  ecs_int_t    val_ref;

  size_t       ielt_ref;
  size_t       ielt_transf;
  size_t       ipos_ref;
  size_t       ipos_transf;

  size_t       elt_nbr_val_max;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (this_table_att == NULL)
    return;

  ecs_table__regle_en_pos(this_table_att);

  nbr_elt_ref = this_table_att->nbr;

  elt_nbr_val_max = 0;

  /* Détermination du nombre maximal de propriétés pour un élément */

  for (ielt_ref = 0; ielt_ref < nbr_elt_ref; ielt_ref++)
    elt_nbr_val_max = ECS_MAX(elt_nbr_val_max,
                              (  this_table_att->pos[ielt_ref + 1]
                               - this_table_att->pos[ielt_ref]));

  /* On suppose qu'une fusion ne concerne qu'au maximum 2 éléments       */
  /* Le nombre maximal de propriétés pour un élément fusionne est donc : */

  elt_nbr_val_max *= 2;

  ECS_MALLOC(pos_tab_transf, nbr_elt_new, ecs_size_t);
  ECS_MALLOC(val_tab_transf, nbr_elt_new, ecs_int_t *);
  ECS_MALLOC(val_tab_unidim, nbr_elt_new * elt_nbr_val_max, ecs_int_t);

  for (ielt_transf = 0; ielt_transf < nbr_elt_new; ielt_transf++) {
    pos_tab_transf[ielt_transf]= 0;
    val_tab_transf[ielt_transf]= &val_tab_unidim[ielt_transf * elt_nbr_val_max];
  }

  nbr_val_transf = 0;

  for (ielt_ref = 0; ielt_ref < nbr_elt_ref; ielt_ref++) {

    num_elt_transf  = vect_transf.val[ielt_ref];

    pos_inf = this_table_att->pos[ielt_ref]     - 1;
    pos_sup = this_table_att->pos[ielt_ref + 1] - 1;

    for (ipos_ref = pos_inf; ipos_ref < pos_sup; ipos_ref++) {

      val_ref = this_table_att->val[ipos_ref];

      ipos_transf = 0;
      while (ipos_transf < pos_tab_transf[num_elt_transf]           &&
             val_tab_transf[num_elt_transf][ipos_transf] != val_ref    )
        ipos_transf++;

      if (ipos_transf == pos_tab_transf[num_elt_transf]) {

        /* C'est une nouvelle valeur : on la stocke */

        val_tab_transf[num_elt_transf][ipos_transf] = val_ref;
        pos_tab_transf[num_elt_transf]++;
        nbr_val_transf++;

      }
      /* else : rien à faire (la valeur a déjà été stockée) */

    } /* Fin : boucle sur les valeurs de l'élément de référence */

  } /* Fin : boucle sur les éléments de référence */

  this_table_att->nbr = nbr_elt_new;

  ECS_REALLOC(this_table_att->pos, nbr_elt_new + 1, ecs_size_t);

  ECS_REALLOC(this_table_att->val, nbr_val_transf, ecs_int_t);

  this_table_att->pos[0] = 1;
  cpt_val_transf = 0;

  for (ielt_transf = 0; ielt_transf < nbr_elt_new; ielt_transf++) {

    this_table_att->pos[ielt_transf + 1]
      = this_table_att->pos[ielt_transf] + pos_tab_transf[ielt_transf];

    for (ipos_transf = 0;
         ipos_transf < pos_tab_transf[ielt_transf];
         ipos_transf++) {

      this_table_att->val[cpt_val_transf++]
        = val_tab_transf[ielt_transf][ipos_transf];
    }
  }

  ECS_FREE(pos_tab_transf);
  ECS_FREE(val_tab_transf);
  ECS_FREE(val_tab_unidim);

  ecs_table__pos_en_regle(this_table_att);
}

/*----------------------------------------------------------------------------
 *  Fonction qui construit les familles à partir
 *   de la liste chaînée de tous les tables de type "attribut"
 *   pour toutes des entités ; ces familles sont ajoutées à la liste
 *   chaînée fournie en argument ;
 *
 *  Elle remplace la table attribut par une table "famille" par entité
 *
 *  Elle détermine aussi :
 *   - le nombre de familles
 *----------------------------------------------------------------------------*/

int *
ecs_table_att__construit_fam(ecs_table_t     **table_att,
                             ecs_famille_t   **vect_fam_tete,
                             int               num_fam_deb,
                             int              *nbr_fam)
{
  size_t         n_elts;
  ecs_int_t      ifam;
  ecs_int_t     *nbr_descr_fam;
  ecs_int_t   **def_fam_descr;

  ecs_tab_int_t  tab_idem;
  ecs_tab_int_t  tab_renum_elt;

  int  *elt_fam = NULL;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (table_att == NULL)
    return NULL;

  if (*table_att == NULL)
    return NULL;

  n_elts = ecs_table__ret_elt_nbr(*table_att);

  /* Tri des valeurs de la table unifiee */
  /*-------------------------------------*/

  _table_att__trie_descr(*table_att);

  ecs_table__regle_en_pos(*table_att);

  tab_renum_elt = _table_att__trie_val(table_att);

  /* Tous les éléments ayant des listes de numéros de descripteurs d'attribut */
  /*  identiques appartiennent à la même famille                              */
  /*  (famille définie par la liste de descripteurs d'attribut)               */

  tab_idem = _table_att__compare_val_pos(*table_att);

  elt_fam = _table__attribue_fam(*table_att,
                                 tab_idem,
                                 &def_fam_descr,
                                 &nbr_descr_fam,
                                 num_fam_deb,
                                 nbr_fam);

  ECS_FREE(tab_idem.val);

  ecs_table__pos_en_regle(*table_att);

  /* Retour à la numérotation initiale des éléments pour les familles */

  _ecs_table__transforme_fam(n_elts,
                             elt_fam,
                             tab_renum_elt);

  ECS_FREE(tab_renum_elt.val);

  /*-------------------------------------------------------------------------*/
  /* Construction des familles à partir :                                    */
  /*  - du vecteur indexe donnant pour chaque numéro de famille              */
  /*     la liste des numéros d'attributs                                    */
  /*  - de la liste chaînée de l'ensemble des descripteurs d'attribut        */
  /*-------------------------------------------------------------------------*/

  *vect_fam_tete = ecs_famille_chaine__cree(def_fam_descr,
                                            nbr_descr_fam,
                                            num_fam_deb,
                                            *nbr_fam,
                                            (*table_att)->descr);

  ecs_table__detruit(table_att);

  for (ifam = 0; ifam < *nbr_fam; ifam++)
    ECS_FREE(def_fam_descr[ifam]);
  ECS_FREE(def_fam_descr);
  ECS_FREE(nbr_descr_fam);

  return elt_fam;
}


/*----------------------------------------------------------------------------
 *  Fonction qui crée la table "groupe" à partir
 *   de la table "famille" et de la liste chaînée des familles
 *----------------------------------------------------------------------------*/

ecs_table_t *
ecs_table_att__cree_att_fam(size_t           n_elts,
                            int             *elt_fam,
                            ecs_famille_t   *famille)
{
  size_t        ifam;
  int           num_fam_max;
  int           nbr_max_att_fam;
  int           nbr_groupe;

  ecs_descr_t    *descr_tete_groupe;

  ecs_tab_int_t   tab_fam;
  ecs_tab_int_t  *tab_groupe_fam;

  ecs_table_t    *table_fam_groupe;

  ecs_table_t    *table_groupe = NULL;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (elt_fam == NULL || famille == NULL)
    return NULL;

  /* Création de la liste des numéros de famille référencés */
  /*  dans le table "famille"                               */
  /*--------------------------------------------------------*/

  tab_fam = _table_att__ret_reference(n_elts, elt_fam);

  num_fam_max = 0;
  for (ifam = 0; ifam < tab_fam.nbr; ifam++) {
    if (tab_fam.val[ifam] > num_fam_max)
      num_fam_max = tab_fam.val[ifam];
  }

  ECS_MALLOC(tab_groupe_fam , num_fam_max, ecs_tab_int_t);

  for (ifam = 0; ifam < (size_t)num_fam_max; ifam++) {
    tab_groupe_fam [ifam].val = NULL;
    tab_groupe_fam [ifam].nbr = 0;
  }

  /* Création d'une liste chaînée de descripteurs pour le table "groupe" */
  /*---------------------------------------------------------------------*/

  ecs_famille_chaine__cree_descr(famille,
                                 tab_fam,
                                 &descr_tete_groupe,
                                 tab_groupe_fam,
                                 &nbr_max_att_fam);

  if (tab_fam.nbr != 0)
    ECS_FREE(tab_fam.val);

  nbr_groupe = 0;
  for (ifam = 0; ifam < (size_t)num_fam_max; ifam++)
    nbr_groupe += tab_groupe_fam[ifam].nbr;

  if (nbr_groupe != 0) {

    /* Création da la table "groupe" */
    /*-------------------------------*/

    table_fam_groupe  = _table_att__transforme_bi_tab(tab_groupe_fam,
                                                      num_fam_max,
                                                      nbr_groupe);
    ecs_table__regle_en_pos(table_fam_groupe);

    table_groupe = _table_att__famille_en_groupe(n_elts,
                                                 elt_fam,
                                                 table_fam_groupe);

    ecs_table__detruit(&table_fam_groupe);

    table_groupe->descr = descr_tete_groupe;

    ecs_table__pos_en_regle(table_groupe);

    for (ifam = 0; ifam < (size_t)num_fam_max; ifam++)
      if (tab_groupe_fam[ifam].val != NULL)
        ECS_FREE(tab_groupe_fam[ifam].val);
  }
  else {

    table_groupe = NULL;

  }

  ECS_FREE(tab_groupe_fam);

  return table_groupe;
}

/*----------------------------------------------------------------------------*/

