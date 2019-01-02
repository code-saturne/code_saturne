/*============================================================================
 *  Définitions des fonctions de base
 *  associées à la structure `ecs_maillage_t' décrivant un maillage.
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

#include "cs_config.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie standard C
 *----------------------------------------------------------------------------*/

#include <assert.h>


/*----------------------------------------------------------------------------
 *  Fichiers `include' système
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/

#include "ecs_def.h"
#include "ecs_fic.h"
#include "ecs_file.h"
#include "ecs_mem.h"
#include "ecs_tab.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage global "Post-Traitement"
 *----------------------------------------------------------------------------*/

#include "ecs_post.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles des paquetages visibles
 *----------------------------------------------------------------------------*/

#include "ecs_famille_chaine.h"
#include "ecs_post.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_table.h"
#include "ecs_table_def.h"
#include "ecs_table_att.h"
#include "ecs_table_post.h"
#include "ecs_maillage_post.h"


/*----------------------------------------------------------------------------
 *  Fichier  `include' du  paquetage courant associé au fichier courant
 *----------------------------------------------------------------------------*/

#include "ecs_maillage.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' privés   du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_maillage_priv.h"


/*============================================================================
 *                              Fonctions privées
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction d'impression d'une table avec position réglée en ASCII
 *----------------------------------------------------------------------------*/

static void
_imprime_coords(FILE               *f,
                size_t              nbr,
                const ecs_coord_t   val[],
                size_t              nbr_imp)
{
  /* Variables locales */

  size_t  ient;

  size_t  ind_ent_1 = 0;
  size_t  ind_ent_2;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(val != NULL);

  ind_ent_2 = ECS_MIN(nbr, nbr_imp);

  /* Impression des valeurs */
  /*========================*/

  while (1) {

    for (ient = ind_ent_1; ient < ind_ent_2; ient++)
      fprintf(f, "%24s %12lu %#12.5E %#12.5E %#12.5E\n",
              " ", (unsigned long)(ient + 1), (double)val[3*ient],
              (double)val[3*ient + 1], (double)val[3*ient + 2]);

    if (ind_ent_2 == nbr)
      break;

    ind_ent_1 = ECS_MAX(nbr - nbr_imp, nbr_imp);

    if (ind_ent_1 > ind_ent_2)

      fprintf(f,
              "%77s", "...........  ...........  ...........\n");

    ind_ent_2 = nbr;

  }
}

/*----------------------------------------------------------------------------
 *  Fonction d'impression des familles des éléments
 *----------------------------------------------------------------------------*/

static void
_dump_elt_fam(FILE       *f,
              size_t      n_elts,
              const int   elt_fam[],
              size_t      nbr_imp)
{
  /* Variables locales */

  size_t  i;

  size_t  i1 = 0;
  size_t  i2;

  /* Instructions */

  assert(elt_fam != NULL);

  i2 = ECS_MIN(n_elts, nbr_imp);

  /* Print values */
  /*==============*/

  while (1) {

    for (i = i1; i < i2; i++)
      fprintf(f, "%54s %10lu %10d\n", " ", (unsigned long)(i+1), elt_fam[i]);

    if (i2 == n_elts)
      break;

    i1 = ECS_MAX(n_elts - nbr_imp, nbr_imp);

    if (i1 > i2)
      fprintf(f,
              "%77s", "..........\n");

    i2 = n_elts;

  }
}

/*----------------------------------------------------------------------------
 *  Fonction d'impression des connectivités supplémentaires
 *----------------------------------------------------------------------------*/

static void
_dump_connect_couples(FILE             *f,
                      size_t            n_connect_couples,
                      const ecs_int_t   connect_couples[],
                      size_t            nbr_imp)
{
  /* Variables locales */

  size_t  i;

  size_t  i1 = 0;
  size_t  i2;

  /* Instructions */

  assert(connect_couples != NULL);

  i2 = ECS_MIN(n_connect_couples, nbr_imp);

  /* Print values */
  /*==============*/

  while (1) {

    for (i = i1; i < i2; i++)
      fprintf(f, "%50s %12ld %12ld\n", " ",
              (long)connect_couples[i*2],
              (long)connect_couples[i*2+1]);

    if (i2 == n_connect_couples)
      break;

    i1 = ECS_MAX(n_connect_couples - nbr_imp, nbr_imp);

    if (i1 > i2)
      fprintf(f,
              "%77s", "............\n");

    i2 = n_connect_couples;

  }
}

/*----------------------------------------------------------------------------
 * Fonction qui selectionne tous les elements ou ceux appartenant a une liste
 * Seules sont concernees les entites de type `entmail_sel'
 *----------------------------------------------------------------------------*/

static bool *
_maillage__selectionne_lst(ecs_maillage_t       *maillage,
                           ecs_entmail_t         entmail_sel,
                           const ecs_tab_int_t  *liste_filtre)
{
  size_t  nbr_elt;
  size_t  ielt;

  bool  *elt_select = NULL;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (maillage->table_def[entmail_sel] != NULL) {

    nbr_elt = ecs_table__ret_elt_nbr(maillage->table_def[entmail_sel]);

    ECS_MALLOC(elt_select, nbr_elt, bool);

    for (ielt = 0; ielt < nbr_elt; ielt++)
      elt_select[ielt] = false;

  }

  if (liste_filtre == NULL) {

    nbr_elt = ecs_table__ret_elt_nbr(maillage->table_def[entmail_sel]);

    for (ielt = 0; ielt < nbr_elt; ielt++)
      elt_select[ielt] = true;

  }
  else {

    for (ielt = 0; ielt < liste_filtre->nbr; ielt++)
      elt_select[liste_filtre->val[ielt]] = true;

  }

  return elt_select;
}

/*----------------------------------------------------------------------------
 *  Fonction qui détermine une nouveau table à partir d'une table de référence
 *   en extrayant de ce dernier les éléments sélectionnés
 *   par le tableau de booléens
 *
 *  Cette fonction renvoie le tableau qui définit les anciens éléments
 *   du table de référence en fonction des nouveaux éléments du table renvoyé
 *----------------------------------------------------------------------------*/

static ecs_tab_int_t
_maillage__extrait_coords(ecs_maillage_t         *maillage_new,
                          const ecs_maillage_t   *maillage_ref,
                          const bool              elt_select[])
{
  size_t  cpt_elt_new;
  size_t  cpt_val_new;
  size_t  nbr_elt_ref;
  size_t  nbr_val_ref;
  size_t  pos_ref_inf;
  size_t  pos_ref_sup;

  size_t  ielt_ref;
  size_t  ipos_ref;

  ecs_tab_int_t  tab_old_new;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  nbr_elt_ref = maillage_ref->n_vertices;
  nbr_val_ref = maillage_ref->n_vertices * 3;

  maillage_new->n_vertices = maillage_ref->n_vertices;
  ECS_MALLOC(maillage_new->vertex_coords, nbr_val_ref, ecs_coord_t);

  tab_old_new.nbr = nbr_elt_ref;
  ECS_MALLOC(tab_old_new.val, tab_old_new.nbr, ecs_int_t);

  cpt_elt_new     = 0;
  cpt_val_new     = 0;

  for (ielt_ref = 0; ielt_ref < nbr_elt_ref; ielt_ref++) {

    if (elt_select[ielt_ref] == true) {

      /* L'élément est à extraire */

      pos_ref_inf = 3 *  ielt_ref;
      pos_ref_sup = 3 * (ielt_ref + 1);

      for (ipos_ref = pos_ref_inf; ipos_ref < pos_ref_sup; ipos_ref++)
        maillage_new->vertex_coords[cpt_val_new++]
          = maillage_ref->vertex_coords[ipos_ref];

      tab_old_new.val[ielt_ref] = cpt_elt_new + 1;

      cpt_elt_new++;

    }
    else {

      tab_old_new.val[ielt_ref] = 0;

    }
  }

  maillage_new->n_vertices = cpt_elt_new * 3;
  ECS_REALLOC(maillage_new->vertex_coords,
              maillage_new->n_vertices*3,
              ecs_coord_t);

  return tab_old_new;
}

/*----------------------------------------------------------------------------
 *  Fonction qui definit de nouvelles entites de maillage principales
 *   par extraction d'une partie des elements
 *   d'une entite de maillage principale donnee
 *  Les elements a extraire sont ceux qui ont un booleen a `true'
 *----------------------------------------------------------------------------*/

static ecs_maillage_t *
_maillage__extrait(ecs_maillage_t  *maillage,
                   ecs_entmail_t    entmail,
                   bool             elt_select[])
{
  size_t          isom;
  ecs_tab_int_t   tab_som_old_new;

  bool            *som_select = NULL;

  ecs_maillage_t  *maillage_new;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(maillage != NULL);

  /* Initialisations */
  /*-----------------*/

  tab_som_old_new.nbr = 0;
  tab_som_old_new.val = NULL;

  maillage_new = ecs_maillage__cree_nodal();

  /* Construction de la connectivité du nouveau maillage extrait */
  /*-------------------------------------------------------------*/

  /* Extraction des elements sélectionnés de l'entite de maillage
     (qui sont renumerotés à partir de 1 mais qui sont toujours
     définis par les sous-éléments non renumérotés)
     Creation d'une table donnant pour chaque ancien numero de l'élément
     avant extraction le nouveau numéro de l'élément selectionné et
     renumeroté a partir de 1  */

  maillage_new->table_def[entmail]
    = ecs_table__extrait(maillage->table_def[entmail],
                         elt_select);

  /* Traitement des sommets */
  /*------------------------*/

  /* Construction de la liste de selection des sommets a extraire */

  ECS_MALLOC(som_select, maillage->n_vertices, bool);

  for (isom = 0; isom < maillage->n_vertices; isom++)
    som_select[isom] = false;

  ecs_table_def__cree_masque(som_select,
                             maillage_new->table_def[entmail]);

  /* Extraction des sommets selectionnés
     Création d'une table donnant
     pour chaque ancien numéro du sommet avant extraction
     le nouveau numéro du sommet selectionné et renumeroté à partir de 1 */

  tab_som_old_new = _maillage__extrait_coords(maillage_new,
                                              maillage,
                                              som_select);

  ECS_FREE(som_select);

  /* Remplacement des anciens numéros des sommets
     par les nouveaux numéros (numérotés à partir de 1)
     dans la définition des éléments */

  ecs_table_def__remplace_ref(maillage_new->table_def[entmail],
                              &tab_som_old_new);

  tab_som_old_new.nbr = 0;
  ECS_FREE(tab_som_old_new.val);

  /* On renvoie le maillage extrait */
  /*--------------------------------*/

  return maillage_new;
}

/*----------------------------------------------------------------------------
 * Concaténation de deux ensembles de sommets.
 *----------------------------------------------------------------------------*/

static void
_maillage__concat_vtx(ecs_maillage_t  *maillage,
                      ecs_maillage_t  *maillage_concat)
{
  size_t    n_vertices_ini;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(maillage   != NULL);
  assert(maillage_concat != NULL);

  n_vertices_ini = maillage->n_vertices;

  maillage->n_vertices += maillage_concat->n_vertices;
  ECS_REALLOC(maillage->vertex_coords,
              maillage->n_vertices * 3,
              ecs_coord_t);

  memcpy(maillage->vertex_coords + (n_vertices_ini*3),
         maillage_concat->vertex_coords,
         maillage_concat->n_vertices * 3 * sizeof(ecs_coord_t));

  /* Remove corresponding data from appended mesh */

  maillage_concat->n_vertices = 0;
  ECS_FREE(maillage_concat->vertex_coords);
}

/*----------------------------------------------------------------------------
 * Concaténation de deux ensembles de connectivités.
 *----------------------------------------------------------------------------*/

static void
_maillage__concat_connect(ecs_maillage_t  *maillage,
                          ecs_maillage_t  *maillage_concat,
                          ecs_entmail_t    entmail)
{
  size_t  n_couples_concat;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(maillage   != NULL);
  assert(maillage_concat != NULL);

  n_couples_concat = maillage_concat->n_connect_couples[entmail];

  if (n_couples_concat > 0) {

    size_t i;
    size_t n_couples = maillage->n_connect_couples[entmail];
    size_t concat_shift = ecs_table__ret_elt_nbr(maillage->table_def[entmail]);
    ecs_int_t  *dest = NULL;
    const ecs_int_t  *src = maillage_concat->connect_couples[entmail];

    /* Update receiving mesh */

    maillage->n_connect_couples[entmail] += n_couples_concat;
    ECS_REALLOC(maillage->connect_couples[entmail],
                maillage->n_connect_couples[entmail]*2,
                ecs_int_t);

    dest = maillage->connect_couples[entmail] + (n_couples*2);

    for (i = 0; i < n_couples_concat*2; i++)
      dest[i] = src[i] + concat_shift;

    /* Remove corresponding data from appended mesh */

    maillage_concat->n_connect_couples[entmail] = 0;
    ECS_FREE(maillage_concat->connect_couples[entmail]);
  }
}

/*----------------------------------------------------------------------------
 *  Fonction réalisant transformation des familles
 *   en fusionnant les familles des éléments qui sont identiquement
 *   transformés par le vecteur de transformation donné.
 *
 *  Cette fonction inclut le prolongement du tableau des familles
 *   d'une liste d'éléments initiale vers une liste à compacter,
 *   ainsi que le compactage en question. Ceci permet de rendre ce
 *   prolongement implicite, et réduire la taille de tableau intermédiaire.
 *
 *  Par construction préalable, on fait l'hypothèse qu'un seul un de chaque
 *   ensemble d'éléments fusionnés porte une famille.
 *----------------------------------------------------------------------------*/

static void
_maillage_elt_fam_fusionne(int                  **elt_fam,
                           size_t                 nbr_elt_old,
                           size_t                 nbr_elt_new,
                           const ecs_tab_int_t    vect_transf)
{
  ecs_int_t    ind_elt_transf;
  int          val_ref;

  size_t       ielt;
  size_t       ielt_ref;

  int         *elt_fam_new = NULL;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (elt_fam == NULL)
    return;

  if (*elt_fam == NULL)
    return;

  ECS_MALLOC(elt_fam_new, nbr_elt_new, int);
  for (ielt = 0; ielt < nbr_elt_new; ielt++)
    elt_fam_new[ielt] = 0;

  for (ielt_ref = 0; ielt_ref < nbr_elt_old; ielt_ref++) {

    ind_elt_transf = vect_transf.val[ielt_ref];
    assert(ind_elt_transf > -1);

    val_ref = (*elt_fam)[ielt_ref];

    if (val_ref != 0 && elt_fam_new[ind_elt_transf] == 0)
      elt_fam_new[ind_elt_transf] = val_ref;

  }

  ECS_FREE(*elt_fam);
  *elt_fam = elt_fam_new;
}

/*----------------------------------------------------------------------------
 *  Update element family numbers when elements are renumbered. Each
 *   initial element has at most one matching element.
 *----------------------------------------------------------------------------*/

static void
_maillage_elt_fam_compacte(int            **elt_fam,
                           ecs_tab_int_t   *tab_old_new)
{
  size_t        ielt;
  ecs_int_t     num_elt_new;

  ecs_int_t     cpt_elt = 0;
  int          *_elt_fam = NULL;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (elt_fam == NULL)
    return;

  if (*elt_fam == NULL)
    return;

  /* Compact numbering */

  _elt_fam = *elt_fam;

  for (ielt = 0; ielt < tab_old_new->nbr; ielt++) {

    num_elt_new = tab_old_new->val[ielt];

    if (num_elt_new != 0) {
      assert(num_elt_new == cpt_elt + 1);
      _elt_fam[num_elt_new - 1] = _elt_fam[ielt];
      cpt_elt += 1;
    }

  }

  /* Update definitions */

  ECS_REALLOC(_elt_fam, cpt_elt, int);

  *elt_fam = _elt_fam;
}

/*----------------------------------------------------------------------------
 *  Suppression des éléments dégénérés
 *----------------------------------------------------------------------------*/

static void
_maillage__nettoie_descend(ecs_maillage_t  *maillage)
{
  ecs_tab_int_t  tab_fac_old_new;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(maillage != NULL);

  if (maillage->table_def[ECS_ENTMAIL_FAC] == NULL)
    return;

  /* Faces */
  /*-------*/

  tab_fac_old_new
    = ecs_table_def__nettoie_fac(maillage->table_def[ECS_ENTMAIL_FAC]);

  if (tab_fac_old_new.nbr != 0) {

    /* Inherit "family" fields */

    assert(maillage->table_att[ECS_ENTMAIL_FAC] == NULL);

    _maillage_elt_fam_compacte(&(maillage->elt_fam[ECS_ENTMAIL_FAC]),
                               &tab_fac_old_new);

    /* Replace references in face definitions */

    if (maillage->table_def[ECS_ENTMAIL_CEL] != NULL)
      ecs_table_def__remplace_ref(maillage->table_def[ECS_ENTMAIL_CEL],
                                  &tab_fac_old_new);

    tab_fac_old_new.nbr = 0;
    ECS_FREE(tab_fac_old_new.val);
  }
}

/*----------------------------------------------------------------------------
 *  Fonction qui construit la liste des faces avec erreur de connectivité
 *  (i.e. qui appartiennent à 2 cellules ou plus vu d'un même côté, ou qui
 *  sont à la fois entrante et sortante pour une cellule).
 *
 *  Un tableau indiquant le type associé à chaque face (0 pour face isolée,
 *  1 ou 2 pour face de bord, 3 pour face interne, et 4 pour tous les autres
 *  cas (faces voyant au moins deux cellules sur un même côté, d'ou erreur
 *  de connectivité) doit être fourni en entrée.
 *----------------------------------------------------------------------------*/

static void
_maillage__liste_fac_err(const ecs_tab_int_t  *typ_fac,
                         ecs_tab_int_t        *liste_fac_erreur)
{
  size_t  cpt_fac_erreur;
  size_t  ifac;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(typ_fac != NULL);
  assert(liste_fac_erreur != NULL);

  /* Initialisations */

  cpt_fac_erreur  = 0;

  /* Première boucle sur les faces : comptage */
  /*------------------------------------------*/

  for (ifac = 0; ifac < typ_fac->nbr; ifac++) {
    if (typ_fac->val[ifac] >= 4)
      cpt_fac_erreur++;
  }

  /* Initialisation et allocation de la liste */

  liste_fac_erreur->nbr = cpt_fac_erreur;
  ECS_MALLOC(liste_fac_erreur->val, liste_fac_erreur->nbr, ecs_int_t);

  /* Seconde boucle sur les faces : remplissage des listes */
  /*-------------------------------------------------------*/

  cpt_fac_erreur  = 0;

  for (ifac = 0; ifac < typ_fac->nbr; ifac++) {
    if (typ_fac->val[ifac] >= 4)
      liste_fac_erreur->val[cpt_fac_erreur++] = ifac;
  }
}

/*----------------------------------------------------------------------------
 *  Fonction qui compte le nombre de faces internes et de bord, ainsi que
 *  le nombre de faces avec erreur de connectivité (i.e. qui appartiennent
 *  à 2 cellules ou plus vu d'un même côté, ou qui sont à la fois entrante
 *  et sortante pour une cellule) et de faces isolées.
 *
 *  Un tableau indiquant le type associé à chaque face (0 pour face isolée,
 *  1 ou 2 pour face de bord, 3 pour face interne, et >= 4 pour tous les autres
 *  cas (faces voyant au moins deux cellules sur un même côté, d'ou erreur
 *  de connectivité) doit être fourni en entrée.
 *----------------------------------------------------------------------------*/

static void
_maillage__compte_typ_fac(const ecs_tab_int_t  *typ_fac,
                          size_t               *nbr_fac_erreur,
                          size_t               *nbr_fac_interne,
                          size_t               *nbr_fac_de_bord,
                          size_t               *nbr_fac_isolee)
{
  size_t  ifac;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(typ_fac != NULL);

  /* Initialisations */

  *nbr_fac_erreur  = 0;
  *nbr_fac_interne = 0;
  *nbr_fac_de_bord = 0;
  *nbr_fac_isolee  = 0;

  /* Boucle sur les faces : comptage */
  /*---------------------------------*/

  for (ifac = 0; ifac < typ_fac->nbr; ifac++) {

    switch(typ_fac->val[ifac]) {

    case 0:
      *nbr_fac_isolee += 1;
      break;

    case 1:
    case 2:
      *nbr_fac_de_bord += 1;
      break;

    case 3:
      *nbr_fac_interne += 1;
      break;

    default:
      *nbr_fac_erreur += 1;

    }
  }

  if (*nbr_fac_erreur != 0) {
    ecs_warn();
    printf(_("There are %lu faces of which one same side belongs\n"
             "to at least 2 cells --> bad connectivity."),
           (unsigned long)(*nbr_fac_erreur));
  }

  if (*nbr_fac_isolee != 0) {
    ecs_warn();
    printf(_("There is/are %lu isolated face(s)\n"),
           (unsigned long)(*nbr_fac_isolee));
  }
}

/*----------------------------------------------------------------------------
 *  Fonction réalisant la transformation d'un vecteur indexé
 *   en appliquant directement le vecteur de transformation donné
 *   sur les valeurs associées à ses éléments
 *----------------------------------------------------------------------------*/

static void
_maillage_renum_connect(ecs_maillage_t       *maillage,
                        ecs_entmail_t         entmail,
                        const ecs_tab_int_t   vect_transf)
{
  size_t       n_couples, i, j;
  ecs_int_t    c1, c2;
  ecs_int_t   *connect;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(maillage != NULL);

  if (maillage->connect_couples[entmail] == NULL)
    return;

  n_couples = maillage->n_connect_couples[entmail];
  connect = maillage->connect_couples[entmail];

  for (i = 0, j = 0; i < n_couples; i++) {
    assert(connect[i*2] > 0 && connect[i*2 + 1] > 0);
    c1 = ECS_ABS(vect_transf.val[connect[i*2] - 1]);
    c2 = ECS_ABS(vect_transf.val[connect[i*2 + 1] - 1]);
    if (c1 != c2 && c1 > 0 && c2 > 0) {
      connect[j++] = c1;
      connect[j++] = c2;
    }
  }

  if (j < n_couples*2) {
    maillage->n_connect_couples[entmail] = j/2;
    ECS_REALLOC(maillage->connect_couples[entmail], j, ecs_int_t);
  }
}

/*============================================================================
 *                             Fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define a new empty mesh structure with nodal connectivity.
 *----------------------------------------------------------------------------*/

ecs_maillage_t *
ecs_maillage__cree_nodal(void)
{
  int ient;
  ecs_maillage_t  *maillage = NULL;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Structure globale du maillage */
  /*-------------------------------*/

  /* Allocation de la structure globale du maillage */

  ECS_MALLOC(maillage, 1, ecs_maillage_t);

  /* Initialisation du type de connectivité par défaut */

  maillage->typ_connect = ECS_MAILLAGE_CONNECT_NODALE;

  /* Entités de maillage */
  /*---------------------*/

  maillage->vertex_coords = NULL;

  for (ient = 0; ient < 2; ient++) {
    maillage->table_def[ient] = NULL;
    maillage->table_att[ient] = NULL;
    maillage->elt_fam[ient] = NULL;
    maillage->n_connect_couples[ient] = 0;
    maillage->connect_couples[ient] = NULL;
    maillage->famille[ient] = NULL;
  }

  return maillage;
}

/*----------------------------------------------------------------------------
 * Free a mesh structure.
 *----------------------------------------------------------------------------*/

void
ecs_maillage__detruit(ecs_maillage_t  **maillage)
{
  int  ient;
  ecs_maillage_t *m = *maillage;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Structures des entités de maillage */
  /*====================================*/

  if (m->vertex_coords != NULL)
    ECS_FREE(m->vertex_coords);

  for (ient = 0; ient < 2; ient++) {
    if (m->table_def[ient] != NULL)
      ecs_table__detruit(&m->table_def[ient]);
    if (m->table_att[ient] != NULL)
      ecs_table__detruit(&m->table_att[ient]);
    if (m->elt_fam[ient] != NULL)
      ECS_FREE(m->elt_fam[ient]);
    if (m->connect_couples[ient] != NULL)
      ECS_FREE(m->connect_couples[ient]);
    if (m->famille[ient] != NULL)
      ecs_famille_chaine__detruit(&(m->famille[ient]));
  }

  ECS_FREE(*maillage);
}

/*----------------------------------------------------------------------------
 *  Fonction imprimant le contenu d'une structure `ecs_maillage_t' donnée
 *   dans le fichier preprocessor_dump.txt
 *----------------------------------------------------------------------------*/

void
ecs_maillage__imprime(const ecs_maillage_t  *maillage,
                      const ecs_int_t        nbr_imp)
{
  const char *nom_typ_c[2] = {
    "ECS_MAILLAGE_CONNECT_NODALE",
    "ECS_MAILLAGE_CONNECT_DESCENDANTE"
  };

  const char *nom_table_def[2] = {
    "FACE_DEFS",
    "CELL_DEFS"
  };

  const char *nom_table_att[2] = {
    "FACE_GROUPS",
    "CELL_GROUPS"
  };

  const char *nom_elt_fam[2] = {
    "FACE_FAMILIES",
    "CELL_FAMILIES"
  };

  const char *nom_connect[2] = {
    "FACE_CONNECTIVITY",
    "CELL_CONNECTIVITY"
  };

  const char *nom_fam[2] = {
    "FACE_FAMILIES",
    "CELL_FAMILIES"
  };

  FILE        *f;
  ecs_int_t    imp_col;
  ecs_int_t    ient;

#define ECS_FCT_IMP_MAILLAGE_FAMILLE       "famille"

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(maillage != NULL);

  imp_col = 0;

  /* Ouverture du fichier d'impression */
  /*===================================*/

  f = fopen("preprocessor_dump.txt", "w");

  /* Message sur la sortie standard */
  /*================================*/

  printf("\n\nMesh structure after reading input\n\n");

  /* Connectivity type */
  /*-------------------*/

  ecs_fic__imprime_val(f, imp_col, "typ_connect",
                       ECS_TYPE_char, nom_typ_c[maillage->typ_connect]);

  /* Vertices */
  /*----------*/

  ecs_fic__imprime_val(f, imp_col, "n_vertices",
                       ECS_TYPE_size_t,
                       &(maillage->n_vertices));

  ecs_fic__imprime_ptr(f, imp_col,
                       "VERTEX_COORDS",
                       maillage->vertex_coords);

  if (maillage->vertex_coords != NULL)
    _imprime_coords(f,
                    maillage->n_vertices,
                    maillage->vertex_coords,
                    nbr_imp);

  /* Element definitions and attributes */
  /*------------------------------------*/

  for (ient = 0; ient < ECS_N_ENTMAIL; ient++) {

    /* Main definitions */

    ecs_fic__imprime_ptr(f, imp_col,
                         nom_table_def[ient],
                         maillage->table_def[ient]);

    if (maillage->table_def[ient] != NULL)
      ecs_table__imprime(maillage->table_def[ient],
                         imp_col+1,
                         nbr_imp,
                         f);

    /* Groups */

    ecs_fic__imprime_ptr(f, imp_col,
                         nom_table_att[ient],
                         maillage->table_att[ient]);

    if (maillage->table_att[ient] != NULL)
      ecs_table__imprime(maillage->table_att[ient],
                         imp_col+1,
                         nbr_imp,
                         f);

    /* Entity families */

    ecs_fic__imprime_ptr(f, imp_col,
                         nom_elt_fam[ient],
                         maillage->elt_fam[ient]);

    if (maillage->elt_fam[ient] != NULL)
      _dump_elt_fam(f,
                    ecs_table__ret_elt_nbr(maillage->table_def[ient]),
                    maillage->elt_fam[ient],
                    nbr_imp);

    /* Additional connectivity */

    ecs_fic__imprime_val(f, imp_col, "n_connect_couples",
                         ECS_TYPE_size_t,
                         &(maillage->n_connect_couples[ient]));

    ecs_fic__imprime_ptr(f, imp_col,
                         nom_connect[ient],
                         maillage->connect_couples[ient]);

    if (maillage->n_connect_couples[ient] != 0)
      _dump_connect_couples(f,
                            maillage->n_connect_couples[ient],
                            maillage->connect_couples[ient],
                            nbr_imp);

  }

  /* Family definitions */
  /*--------------------*/

  for (ient = 0; ient < ECS_N_ENTMAIL; ient++) {

    /* Impression du pointeur sur une entité principale */

    ecs_fic__imprime_ptr(f, imp_col,
                         nom_fam[ient],
                         (void *)maillage->famille[ient]);

    if (maillage->famille[ient] != NULL)
      ecs_famille_chaine__imprime(maillage->famille[ient],
                                  imp_col + 1,
                                  f);
  }

  /* Close dump file */
  /*-----------------*/

  fclose(f);
}

/*----------------------------------------------------------------------------
 *  Fonction qui retourne le type d'entité de plus grande dimension
 *   contenue dans une structure `ecs_maillage_t'
 *----------------------------------------------------------------------------*/

ecs_entmail_t
ecs_maillage__ret_entmail_max(const ecs_maillage_t  *maillage)
{
  ecs_entmail_t entmail_max = ECS_ENTMAIL_NONE;
  ecs_int_t ient;

  for (ient = 0; ient < ECS_N_ENTMAIL; ient++) {
    if (maillage->table_def[ient] != NULL) {
      if (ecs_table__ret_elt_nbr(maillage->table_def[ient]) > 0)
        entmail_max = (ecs_entmail_t) ient;
    }
  }

  return entmail_max;
}

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie la taille en octets d'une structure `ecs_maillage_t'.
 *----------------------------------------------------------------------------*/

float
ecs_maillage__ret_taille(const ecs_maillage_t  *maillage)
{
  size_t       taille;
  ecs_int_t    ient;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(maillage != NULL);

  taille = sizeof(*maillage);

  for (ient = 0; ient < ECS_N_ENTMAIL; ient++) {

    if (maillage->table_def[ient] != NULL)
      taille += ecs_table__ret_taille(maillage->table_def[ient]);

    if (maillage->table_att[ient] != NULL)
      taille += ecs_table__ret_taille(maillage->table_att[ient]);

    if (maillage->elt_fam[ient] != NULL)
      taille += (  ecs_table__ret_elt_nbr(maillage->table_def[ient])
                 * sizeof(int));

    if (maillage->n_connect_couples[ient] != 0)
      taille += maillage->n_connect_couples[ient] * sizeof(ecs_int_t);

    if (maillage->famille[ient] != NULL)
      taille += ecs_famille_chaine__ret_taille(maillage->famille[ient]);
  }

  return (float)taille;
}

/*----------------------------------------------------------------------------
 *  Suppression des sommets ne participant pas à la connectivité
 *   et fusion des éléments surfaciques confondus éventuels
 *----------------------------------------------------------------------------*/

void
ecs_maillage__nettoie_nodal(ecs_maillage_t  *maillage)
{
  size_t         nbr_elt_new;
  ecs_tab_int_t  vect_transf;
  ecs_tab_int_t  signe_elt;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Suppression des sommets inutiles */

  if (maillage->vertex_coords != NULL)
    ecs_table_def__nettoie_nodal(&(maillage->n_vertices),
                                 &(maillage->vertex_coords),
                                 maillage->table_def[ECS_ENTMAIL_FAC],
                                 maillage->table_def[ECS_ENTMAIL_CEL]);

  /*
    Fusion d'éléments surfaciques confondus éventuels (issus par exemple,
    de l'utilisation conjointe de references de faces et faces de bord
    sous Simail, double `surface coating' sous I-DEAS, ...)
  */

  /*--------------------------------------------------------------------------*/
  /* Tri et compactage des definitions des elements                           */
  /*                                                                          */
  /* Determination du vecteur de transformation permettant de passer          */
  /*  de la liste initiale des elements                                       */
  /*  a une liste ordonnee et compactee des elements                          */
  /*--------------------------------------------------------------------------*/

  nbr_elt_new = 0;

  signe_elt.nbr = 0;
  signe_elt.val = NULL;

  if (maillage->table_def[ECS_ENTMAIL_FAC] == NULL)
    return;

  vect_transf = ecs_table_def__fusionne(maillage->table_def[ECS_ENTMAIL_FAC],
                                        &nbr_elt_new,
                                        &signe_elt);

  /* Tables de type attribut ou famille */
  /*------------------------------------*/

  ecs_table_att__fusionne(maillage->table_att[ECS_ENTMAIL_FAC],
                          nbr_elt_new,
                          vect_transf);

  assert(maillage->elt_fam[ECS_ENTMAIL_FAC] == NULL);

  /* Table de type connectivité supplémentaire. */
  /*--------------------------------------------*/

  if (maillage->n_connect_couples[ECS_ENTMAIL_FAC] != 0)
    _maillage_renum_connect(maillage,
                            ECS_ENTMAIL_FAC,
                            vect_transf);

  ECS_FREE(signe_elt.val);
  ECS_FREE(vect_transf.val);
}

/*----------------------------------------------------------------------------
 *  Correction si nécessaire de l'orientation des éléments en
 *   connectivité nodale.
 *
 *  La liste de cellules avec erreur est optionnelle.
 *----------------------------------------------------------------------------*/

void
ecs_maillage__orient_nodal(ecs_maillage_t    *maillage,
                           ecs_tab_int_t     *liste_cel_err,
                           bool               correc_orient)
{
  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (maillage->vertex_coords == NULL)
    return;

  ecs_table_def__orient_nodal(maillage->vertex_coords,
                              maillage->table_def[ECS_ENTMAIL_FAC],
                              maillage->table_def[ECS_ENTMAIL_CEL],
                              liste_cel_err,
                              correc_orient);
}

/*----------------------------------------------------------------------------
 *  Fonction qui assigne la tête de la liste chaînée des familles donnée
 *   à la structure de maillage donnée
 *----------------------------------------------------------------------------*/

void
ecs_maillage__definit_famille(ecs_maillage_t   *maillage,
                              ecs_famille_t    *vect_famille[2])
{
  int ient;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(maillage != NULL);

  for (ient = 0; ient < ECS_N_ENTMAIL; ient++)
    maillage->famille[ient] = vect_famille[ient];
}

/*----------------------------------------------------------------------------
 *  Fonction réalisant, à partir d'une connectivité de maillage donnée,
 *   la connectivité descendante du maillage
 *----------------------------------------------------------------------------*/

void
ecs_maillage__connect_descend(ecs_maillage_t * maillage)
{
  size_t         nbr_elt_new;
  ecs_tab_int_t  vect_transf;
  ecs_tab_int_t  signe_elt;

  ecs_size_t nbr_fac_old = 0;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(maillage != NULL);
  assert(maillage->typ_connect == ECS_MAILLAGE_CONNECT_NODALE);

  maillage->typ_connect = ECS_MAILLAGE_CONNECT_DESCENDANTE;

  if (maillage->table_def[ECS_ENTMAIL_CEL] == NULL)
    return;

  /* Decompose cells into faces */
  /*----------------------------*/

  if (maillage->table_def[ECS_ENTMAIL_FAC] != NULL)
    nbr_fac_old = ecs_table__ret_elt_nbr(maillage->table_def[ECS_ENTMAIL_FAC]);

  ecs_table_def__decompose_cel(&maillage->table_def[ECS_ENTMAIL_FAC],
                               maillage->table_def[ECS_ENTMAIL_CEL]);

  assert(maillage->table_att[ECS_ENTMAIL_FAC] == NULL);

  /* Merge coincident vertices (update face connectivity ) */
  /*-------------------------------------------------------*/

  ecs_table_def__nettoie_som_fac(&(maillage->n_vertices),
                                 &(maillage->vertex_coords),
                                 maillage->table_def[ECS_ENTMAIL_FAC]);

  /* Merge faces with the same definition */
  /*--------------------------------------*/

  nbr_elt_new = 0;

  signe_elt.nbr = 0;
  signe_elt.val = NULL;

  vect_transf = ecs_table_def__fusionne(maillage->table_def[ECS_ENTMAIL_FAC],
                                        &nbr_elt_new,
                                        &signe_elt);

  /* Application du vecteur de transformation sur les autres tables */
  /*----------------------------------------------------------------*/

  assert(maillage->table_att[ECS_ENTMAIL_FAC] == NULL);

  _maillage_elt_fam_fusionne(&(maillage->elt_fam[ECS_ENTMAIL_FAC]),
                             nbr_fac_old,
                             nbr_elt_new,
                             vect_transf);

  _maillage_renum_connect(maillage,
                          ECS_ENTMAIL_FAC,
                          vect_transf);

  /* Application du vecteur de transformation et du signe des elements
     sur la definition des cellules */

  ecs_table__renumerote(maillage->table_def[ECS_ENTMAIL_CEL],
                        vect_transf,
                        signe_elt);

  ECS_FREE(signe_elt.val);
  ECS_FREE(vect_transf.val);

  /* Remove degenerate (empty) faces if present */

  _maillage__nettoie_descend(maillage);
}

/*----------------------------------------------------------------------------
 *  Fonction réalisant le tri des éléments suivant leur type géométrique
 *  La fonction affiche le nombre d'éléments par type géométrique
 *----------------------------------------------------------------------------*/

void
ecs_maillage__trie_typ_geo(ecs_maillage_t  *maillage)
{
  int            ient;
  int            dim_elt[2] = {2, 3};
  ecs_tab_int_t  vect_renum;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(maillage != NULL);

  for (ient = 0; ient < ECS_N_ENTMAIL; ient++) {

    /* Tri des types géometriques des éléments (si nécessaire) */
    /*---------------------------------------------------------*/

    vect_renum.nbr = 0;
    vect_renum.val = NULL;

    if (maillage->table_def[ient] != NULL)

      vect_renum = ecs_table_def__trie_typ(maillage->table_def[ient],
                                           dim_elt[ient]);

    /* Application du vecteur de renumerotation sur les autres tables */
    /*----------------------------------------------------------------*/

    if (vect_renum.val != NULL) {

      /* Inversion du tableau de renumerotation */

      ecs_tab_int__inverse(&vect_renum);

      /* Traitement de la table representant les définitions */

      ecs_table__transforme_pos(maillage->table_def[ient],
                                vect_renum.nbr,
                                vect_renum);

      /* Traitement des tables "attribut" */

      ecs_table__transforme_pos(maillage->table_att[ient],
                                vect_renum.nbr,
                                vect_renum);

      assert(maillage->elt_fam[ient] == NULL);

      /* Traitement des connectivités supplémentaires */

      _maillage_renum_connect(maillage,
                              ient,
                              vect_renum);

      ECS_FREE(vect_renum.val);

    } /* Fin : si le vecteur de renumerotation n'est pas NULL */
  }
}

/*----------------------------------------------------------------------------
 *  Fonction qui définit un nouveau maillage
 *   par extraction d'une partie du maillage donné
 *
 *  Les éléments à extraire doivent être tous de même dimension :
 *  cellules ou faces ou arêtes ou sommets
 *----------------------------------------------------------------------------*/

ecs_maillage_t *
ecs_maillage__extrait(ecs_maillage_t       *maillage,
                      ecs_entmail_t         entmail_sel,
                      const ecs_tab_int_t  *liste_filtre)
{
  bool            *elt_select = NULL;
  ecs_maillage_t  *maillage_new = NULL;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(maillage != NULL);

  assert(   maillage->typ_connect == ECS_MAILLAGE_CONNECT_NODALE
         || entmail_sel == ECS_ENTMAIL_FAC);

  /* Construction de la liste des éléments du maillage d'origine à extraire */
  /*------------------------------------------------------------------------*/

  elt_select = _maillage__selectionne_lst(maillage,
                                          entmail_sel,
                                          liste_filtre);

  /* Extraction des éléments sélectionnés du maillage d'origine */
  /*  pour l'ensemble des entités de maillage                   */
  /*------------------------------------------------------------*/

  maillage_new = _maillage__extrait(maillage,
                                    entmail_sel,
                                    elt_select);

  if (elt_select != NULL)
    ECS_FREE(elt_select);

  return maillage_new;
}

/*----------------------------------------------------------------------------
 *  Fonction qui concatène dans un maillage récepteur donné,
 *   un maillage à concaténer donné.
 *
 *  Le maillage à concaténer est détruit.
 *----------------------------------------------------------------------------*/

void
ecs_maillage__concatene_nodal(ecs_maillage_t  *maillage_recept,
                              ecs_maillage_t  *maillage_concat)
{
  ecs_int_t  ient;
  ecs_int_t  nbr_elt_concat = 0;
  ecs_int_t  nbr_elt_recept = 0;
  ecs_int_t  nbr_som_recept = 0;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(maillage_recept != NULL);
  assert(maillage_concat != NULL);

  assert(maillage_recept->typ_connect == ECS_MAILLAGE_CONNECT_NODALE);
  assert(maillage_concat->typ_connect == ECS_MAILLAGE_CONNECT_NODALE);

  nbr_som_recept = maillage_recept->n_vertices;

  /* Concaténation des sommets */
  /*----------------------------*/

  _maillage__concat_vtx(maillage_recept,
                        maillage_concat);

  /* On doit avoir des attributs et non des familles à ce stade */

  /* Concaténation des éléments */
  /*----------------------------*/

  for (ient = 0; ient < ECS_N_ENTMAIL; ient++) {

    assert(   maillage_recept->elt_fam[ient] == NULL
           && maillage_concat->elt_fam[ient] == NULL);

    /* Décalage des références des sommets */

    if (maillage_concat->table_def[ient] != NULL) {

      ecs_table__incremente_val(maillage_concat->table_def[ient],
                                nbr_som_recept);

      nbr_elt_recept = ecs_table__ret_elt_nbr(maillage_recept->table_def[ient]);
      nbr_elt_concat = ecs_table__ret_elt_nbr(maillage_concat->table_def[ient]);

      if (nbr_elt_recept != 0) {

        /* Définitions */

        ecs_table__concatene(&maillage_recept->table_def[ient],
                             &maillage_concat->table_def[ient],
                             nbr_elt_recept,
                             nbr_elt_concat);

        ecs_table__detruit(&maillage_concat->table_def[ient]);

        /* Groupes */

        ecs_table__concatene(&maillage_recept->table_att[ient],
                             &maillage_concat->table_att[ient],
                             nbr_elt_recept,
                             nbr_elt_concat);

        ecs_table__detruit(&maillage_concat->table_att[ient]);

      }
      else {

        maillage_recept->table_def[ient] = maillage_concat->table_def[ient];
        maillage_recept->table_att[ient] = maillage_concat->table_att[ient];
        maillage_recept->elt_fam[ient] = maillage_concat->elt_fam[ient];

        maillage_concat->table_def[ient] = NULL;
        maillage_concat->table_att[ient] = NULL;
        maillage_concat->elt_fam[ient] = NULL;

      }

      /* Connectivités supplémentaires */

      _maillage__concat_connect(maillage_recept,
                                maillage_concat,
                                ient);

    }  /* else : rien à faire */

  }  /* Fin : boucle sur `ient' */

  ecs_maillage__detruit(&maillage_concat);
}

/*----------------------------------------------------------------------------
 *  Fonction qui construit la liste des cellules attachées à une liste
 *  de faces fournie en argument.
 *----------------------------------------------------------------------------*/

ecs_tab_int_t
ecs_maillage__liste_cel_fac(ecs_maillage_t       *maillage,
                            const ecs_tab_int_t   liste_fac)
{
  size_t  nbr_fac = 0;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(maillage != NULL);
  assert(maillage->table_def[ECS_ENTMAIL_CEL] != NULL);
  assert(maillage->table_def[ECS_ENTMAIL_FAC] != NULL);

  nbr_fac = ecs_table__ret_elt_nbr
             (maillage->table_def[ECS_ENTMAIL_FAC]);

  return ecs_table_def__liste_cel_fac(nbr_fac,
                                      maillage->table_def[ECS_ENTMAIL_CEL],
                                      liste_fac);
}

/*----------------------------------------------------------------------------
 *  Fonction qui calcule les coordonnées min et max du domaine
 *----------------------------------------------------------------------------*/

void
ecs_maillage__calc_coo_ext(ecs_maillage_t  *maillage)
{
  size_t  icoo, ipos, isom, nbr;

  ecs_coord_t  coo_min[3], coo_max[3];

  const ecs_coord_t  *vertex_coords = NULL;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(maillage != NULL);
  assert(maillage->vertex_coords != NULL);

  nbr = maillage->n_vertices;
  vertex_coords = maillage->vertex_coords;

  if (nbr < 1)
    return;

  ipos = 0;

  for (icoo = 0; icoo < 3; icoo++) {
    coo_min[icoo] = vertex_coords[ipos + icoo];
    coo_max[icoo] = vertex_coords[ipos + icoo];
  }

  for (isom = 1; isom < nbr; isom++) {

    ipos = 3 * isom;

    for (icoo = 0; icoo < 3; icoo++) {
      if (vertex_coords[ipos + icoo] < coo_min[icoo])
        coo_min[icoo] = vertex_coords[ipos + icoo];
      else if (vertex_coords[ipos + icoo] > coo_max[icoo])
        coo_max[icoo] = vertex_coords[ipos + icoo];
    }

  }

  printf(_("\n  Domain coordinate extents:\n\n"));

  printf("  [% 10.5e, % 10.5e, % 10.5e]\n",
         coo_min[0], coo_min[1], coo_min[2]);
  printf("  [% 10.5e, % 10.5e, % 10.5e]\n",
         coo_max[0], coo_max[1], coo_max[2]);
}

/*----------------------------------------------------------------------------
 *  Fonction qui transforme les attributs en familles
 *----------------------------------------------------------------------------*/

void
ecs_maillage__cree_famille(ecs_maillage_t  *maillage)
{
  int             ient;
  ecs_int_t       num_fam_deb;
  ecs_int_t       ifam_ent;

  ecs_table_t     *table_att = NULL;

  ecs_famille_t  *vect_famille[2] = {NULL, NULL};
  int             nbr_fam_ent[2] = {0, 0};

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Conctruction des familles */

  num_fam_deb = 1;

  /*------------------------*/
  /* Boucle sur les entites */
  /*------------------------*/

  for (ient = ECS_ENTMAIL_CEL; ient >= ECS_ENTMAIL_FAC; ient--) {

    /* Recuperation de l'adresse de la table des groupes */

    table_att = maillage->table_att[ient];

    if (table_att != NULL) {

      maillage->table_att[ient] = NULL;

      /*------------------------------------*/
      /* Construction des familles à partir */
      /*  des tables "attribut" de l'entité */
      /*------------------------------------*/

      maillage->elt_fam[ient]
        = ecs_table_att__construit_fam(&table_att,
                                       &(vect_famille[ient]),
                                       num_fam_deb,
                                       &(nbr_fam_ent[ient]));

      num_fam_deb += nbr_fam_ent[ient];

    } /* Fin : si il y a des tables "attribut" pour cette entité */

  } /* Fin : boucle sur les entités de maillage */

  /* Récupération des valeurs */

  for (ifam_ent = 0; ifam_ent < ECS_N_ENTMAIL; ifam_ent++)
    maillage->famille[ifam_ent] = vect_famille[ifam_ent];
}

/*----------------------------------------------------------------------------
 *  Fonction qui détruit les familles
 *----------------------------------------------------------------------------*/

void
ecs_maillage__detruit_famille(ecs_maillage_t  *maillage)
{
  int  ient;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Destruction des tables "famille" et des familles */
  /*-------------------------------------------------*/

  for (ient = 0; ient < ECS_N_ENTMAIL; ient++) {

    ECS_FREE(maillage->elt_fam[ient]);

    ecs_famille_chaine__detruit(&maillage->famille[ient]);
  }
}

/*----------------------------------------------------------------------------
 *  Fonction qui construit les attributs "groupe" à partir des familles
 *----------------------------------------------------------------------------*/

void
ecs_maillage__cree_attributs(ecs_maillage_t  *maillage)
{
  int  ient;

  ecs_famille_t  *famille = NULL;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Concaténation des familles des différentes entités */

  famille = NULL;

  for (ient = ECS_ENTMAIL_CEL; ient >= ECS_ENTMAIL_FAC; ient--) {

    if (maillage->famille[ient] != NULL) {

      ecs_famille_chaine__ajoute(&famille, maillage->famille[ient]);

      maillage->famille[ient] = NULL;

    }
  }

  if (famille != NULL) {

    for (ient = ECS_ENTMAIL_CEL; ient >= ECS_ENTMAIL_FAC; ient--) {

      if (maillage->elt_fam[ient] != NULL ) {

        assert(maillage->table_att[ient] == NULL);

        /* Création de la table "groupe" */

        maillage->table_att[ient]
          = ecs_table_att__cree_att_fam
              (ecs_table__ret_elt_nbr(maillage->table_def[ient]),
               maillage->elt_fam[ient],
               famille);

        /* Libération des tableaux de familles */

        ECS_FREE(maillage->elt_fam[ient]);

      }
    }
  }

  ecs_famille_chaine__detruit(&famille);
}

/*----------------------------------------------------------------------------
 *  Fonction qui supprime les attributs "groupe"
 *----------------------------------------------------------------------------*/

void
ecs_maillage__supprime_attributs(ecs_maillage_t  *maillage)
{
  int  ient;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for (ient = 0; ient < ECS_N_ENTMAIL; ient++) {

    if (maillage->table_att[ient] != NULL)
      ecs_table__detruit(&(maillage->table_att[ient]));
  }
}

/*----------------------------------------------------------------------------
 *  Vérification d'un maillage
 *----------------------------------------------------------------------------*/

bool
ecs_maillage__verif(ecs_maillage_t  *maillage,
                    ecs_post_t      *cas_post)
{
  size_t  nbr_cel, nbr_som;
  size_t  nbr_fac_erreur, nbr_fac_interne, nbr_fac_de_bord, nbr_fac_isolee;

  ecs_tab_int_t  typ_fac_cel;
  ecs_tab_int_t  liste_fac_erreur;

  bool  bool_coherent;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(maillage != NULL);

  bool_coherent = true;

  /* On vérifie qu'on est bien en connectivité descendante */

  assert(maillage->typ_connect == ECS_MAILLAGE_CONNECT_DESCENDANTE);

  assert(maillage->vertex_coords != NULL);
  assert(maillage->table_def[ECS_ENTMAIL_FAC] != NULL);

  if (maillage->table_def[ECS_ENTMAIL_CEL] == NULL)
    return false;

  /* Détermination du nombre de cellules et de faces */

  nbr_cel = ecs_table__ret_elt_nbr(maillage->table_def[ECS_ENTMAIL_CEL]);
  nbr_som = maillage->n_vertices;

  /* Determination du type de connectivité associé à chaque face */

  typ_fac_cel
    = ecs_table_def__typ_fac_cel(maillage->table_def[ECS_ENTMAIL_CEL],
                                 maillage->table_def[ECS_ENTMAIL_FAC]);

  _maillage__compte_typ_fac(&typ_fac_cel,
                            &nbr_fac_erreur,
                            &nbr_fac_interne,
                            &nbr_fac_de_bord,
                            &nbr_fac_isolee);

  if (nbr_fac_erreur > 0)
    bool_coherent = false;

  /* Affichage des infos sur les maillage */
  /*--------------------------------------*/

  printf(_("\n\nMain mesh properties\n"
           "--------------------\n\n"));

  printf(_("  Number of cells:                            %10d\n"
           "  Number of internal faces:                   %10d\n"),
         (int)nbr_cel,
         (int)nbr_fac_interne);


  printf(_("  Number of boundary faces:                   %10d\n"),
         (int)nbr_fac_de_bord);


  if (nbr_som != 0)
    printf(_("  Number of vertices:                         %10d\n"),
           (int)nbr_som);

  /* Construction des listes de faces en cas de post-traitement */

  _maillage__liste_fac_err(&typ_fac_cel,
                           &liste_fac_erreur);

  typ_fac_cel.nbr = 0;
  ECS_FREE(typ_fac_cel.val);

  /* Affichage des faces avec connectivité excessive */

  if (liste_fac_erreur.nbr > 0) {

    /* En case de faces avec erreur de connectivité,
       on effectue un post traitement pour analyse de problème */

    ecs_maillage_post__ecr_fac_liste(_("Connectivity Error Faces"),
                                     maillage,
                                     liste_fac_erreur,
                                     ECS_POST_TYPE_ERREUR,
                                     cas_post);

    ECS_FREE(liste_fac_erreur.val);
    liste_fac_erreur.nbr = 0;

  }

  return bool_coherent;
}

/*----------------------------------------------------------------------------*/


