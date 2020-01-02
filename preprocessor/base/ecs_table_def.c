/*============================================================================
 *  Definitions des fonctions
 *   associees a la structure `ecs_table_t' decrivant un table
 *   et propres aux tables principaux de type "definition"
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


/*============================================================================
 *                                 Visibilite
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie standard C
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdlib.h>
#include <string.h>


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/

#include "ecs_elt_typ_liste.h"
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


/*----------------------------------------------------------------------------
 *  Fichier  `include' du  paquetage courant associe au fichier courant
 *----------------------------------------------------------------------------*/

#include "ecs_table_def.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' prives   du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_table_priv.h"


/*============================================================================
 *                       Macros globales au fichier
 *============================================================================*/

/* Longueur maximale du nom d'un type d'élément (+1 pour le `\0' !) */
#define ECS_LOC_LNG_MAX_NOM_TYP    11

#if !defined(FLT_MAX)
#define FLT_MAX HUGE_VAL
#endif

#define ECS_LOC_PRODUIT_VECTORIEL(prod_vect, vect1, vect2) ( \
prod_vect[0] = vect1[1] * vect2[2] - vect2[1] * vect1[2],   \
prod_vect[1] = vect2[0] * vect1[2] - vect1[0] * vect2[2],   \
prod_vect[2] = vect1[0] * vect2[1] - vect2[0] * vect1[1]   )

#define ECS_LOC_PRODUIT_SCALAIRE(vect1, vect2)                        ( \
  vect1[0] * vect2[0] + vect1[1] * vect2[1] + vect1[2] * vect2[2] )

#define ECS_LOC_MODULE(vect) \
     sqrt(vect[0] * vect[0] + vect[1] * vect[1] + vect[2] * vect[2])

#define ECS_LOC_DETERMINANT(vect1, vect2, vect3) ( \
   ((vect1[1] * vect2[2] - vect2[1] * vect1[2]) * vect3[0]) \
 + ((vect2[0] * vect1[2] - vect1[0] * vect2[2]) * vect3[1]) \
 + ((vect1[0] * vect2[1] - vect2[0] * vect1[1]) * vect3[2]) )

/*============================================================================
 *                              Fonctions privees
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction qui met à jour la définition faces -> sommets en cas
 *  de fusion de sommets.
 *----------------------------------------------------------------------------*/

static void
_table_def__maj_fac_som(ecs_table_t          *table_def_fac,
                        const ecs_tab_int_t  *tab_som_old_new)
{
  size_t     cpt_val;
  size_t     nbr_fac;
  size_t     nbr_fac_mod;
  size_t     nbr_val_ini, nbr_val_fin;
  size_t     nbr_som_fac;

  size_t     num_som, num_som_prev;

  size_t     ind_fac;
  size_t     ind_fac_mod;
  size_t     ind_som;

  size_t     ipos_deb, ipos_fin;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Initialisations */
  /* --------------- */

  nbr_fac = table_def_fac->nbr;

  nbr_val_ini = table_def_fac->pos[nbr_fac];

  /* Mise à jour de la définition des faces */
  /*----------------------------------------*/

  for (ind_som = 0;
       ind_som < table_def_fac->pos[table_def_fac->nbr] - 1;
       ind_som++)

    table_def_fac->val[ind_som]
      = tab_som_old_new->val[table_def_fac->val[ind_som] - 1];

  /* Suppression de sommets confondus de la définition des faces */
  /*-------------------------------------------------------------*/

  nbr_fac_mod = 0;

  cpt_val = 0;
  ipos_deb = 1;

  for (ind_fac = 0; ind_fac < nbr_fac; ind_fac++) {

    ind_fac_mod = 0;

    ipos_fin = table_def_fac->pos[ind_fac + 1] - 1;

    nbr_som_fac = ipos_fin - ipos_deb;

    num_som_prev = table_def_fac->val[ipos_deb + nbr_som_fac - 1];

    for (ind_som = 0; ind_som < nbr_som_fac; ind_som++) {

      num_som = table_def_fac->val[ipos_deb + ind_som];

      if (num_som != num_som_prev) {
        num_som_prev = num_som;
        table_def_fac->val[cpt_val++] = num_som;
      }
      else
        ind_fac_mod = 1;
    }

    table_def_fac->pos[ind_fac + 1] = cpt_val + 1;

    ipos_deb = ipos_fin;

    nbr_fac_mod += ind_fac_mod;
  }

  nbr_val_fin = table_def_fac->pos[nbr_fac];

  assert(nbr_val_fin <= nbr_val_ini);

  if (nbr_val_fin != nbr_val_ini) {

    ECS_REALLOC(table_def_fac->val, nbr_val_fin, ecs_int_t);
    printf(_("\nMesh verification:\n\n"
             "  %d faces modified due to merged vertices.\n"),
           (int)nbr_fac_mod);
  }
}

/*----------------------------------------------------------------------------
 *  Fonction qui transforme un tableau d'équivalence en liste chaînée simple
 *----------------------------------------------------------------------------*/

static void
_table_def__transf_equiv(size_t          nbr_som,
                         ecs_tab_int_t  *tab_equiv_som)
{
  size_t     ind_som;
  ecs_int_t  ind_som_min;
  ecs_int_t  ind_som_max;
  ecs_int_t  ind_som_tmp;

  ecs_tab_int_t    tab_equiv_prec;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  tab_equiv_prec.nbr = nbr_som;
  ECS_MALLOC(tab_equiv_prec.val, tab_equiv_prec.nbr, ecs_int_t);

  for (ind_som = 0; ind_som < nbr_som; ind_som++)
    tab_equiv_prec.val[ind_som] = -1;

  for (ind_som = 0; ind_som < nbr_som; ind_som++) {

    if (tab_equiv_som->val[ind_som] != -1) {

      ind_som_min = ind_som;
      ind_som_max = tab_equiv_som->val[ind_som];

      assert (ind_som_min < ind_som_max);

      while (   ind_som_min != ind_som_max
             && tab_equiv_prec.val[ind_som_max] != ind_som_min) {

        /*
          On parcourt la liste inverse correspondant à ind_som_max jusqu'au
          point d'insertion (si pas de point precedent dans la chaine,
          tab_equiv_prec.val[ind_som_max] = -1).
        */

        while (tab_equiv_prec.val[ind_som_max] > ind_som_min)
          ind_som_max = tab_equiv_prec.val[ind_som_max];

        ind_som_tmp = tab_equiv_prec.val[ind_som_max];

        /*
          Si l'on est en début de chaîne, on branche la liste inverse
          correspondant à ind_som_min au début de celle correspondant
          à ind_som_max. Sinon, on doit reboucler.
        */

        tab_equiv_prec.val[ind_som_max] = ind_som_min;

        if (ind_som_tmp != -1) {
          ind_som_max = ind_som_min;
          ind_som_min = ind_som_tmp;
        }

      }
    }
  }

  for (ind_som = 0; ind_som < nbr_som; ind_som++) {

    if (tab_equiv_prec.val[ind_som] != -1)
      tab_equiv_som->val[tab_equiv_prec.val[ind_som]] = ind_som;

  }

  tab_equiv_prec.nbr = 0;
  ECS_FREE(tab_equiv_prec.val);
}

/*----------------------------------------------------------------------------
 *  Fonction qui fusionne des sommets équivalents
 *
 *  Remarque : le tableau d'équivalence (fusion) des sommets est construit de
 *             manière à ce qu'à un sommet ne comportant pas d'équivalent
 *             où de plus petit indice parmi ses équivalents corresponde la
 *             valeur -1, alors qu'un un sommet possédant des équivalents de
 *             plus petit indice corresponde le plus grand indice parmi ces
 *             équivalents (ce qui constitue une sorte de liste chaînée).
 *----------------------------------------------------------------------------*/

static ecs_tab_int_t
_table_def__fusion_som(size_t          *n_vertices,
                       ecs_coord_t    **coords,
                       ecs_tab_int_t   *tab_equiv_som)
{
  size_t  ind_som, nbr_som_old, nbr_som_new, nbr_som_fus;
  ecs_int_t  ind_som_loc, ind_som_tmp, icoo;
  ecs_coord_t  som_poids;
  double  som_coord[3];

  ecs_tab_int_t  tab_som_old_new;

  /* Initialization */
  /* -------------- */

  nbr_som_old = *n_vertices;

  printf(_("\n  Merging vertices:\n"));

  /* Transform vertex equivalences array into simple linked list. */

  _table_def__transf_equiv(nbr_som_old, tab_equiv_som);

  printf(_("    Initial number of vertices        : %10d\n"), (int)nbr_som_old);

  tab_som_old_new.nbr = nbr_som_old;
  ECS_MALLOC(tab_som_old_new.val, tab_som_old_new.nbr, ecs_int_t);

  /* Initialize vertex renumbering array */

  for (ind_som = 0; ind_som < nbr_som_old; ind_som++)
    tab_som_old_new.val[ind_som] = 0;

  /* Main loop on vertices */
  /*-----------------------*/

  nbr_som_new = 0;

  for (ind_som = 0; ind_som < nbr_som_old; ind_som++) {

    /* If the vertex has not been handled yet */

    if (tab_som_old_new.val[ind_som] == 0) {

      /* Initialize vertex */

      som_poids = 0.0;

      for (icoo = 0; icoo < 3; icoo++)
        som_coord[icoo] = 0.0;

      ind_som_loc = ind_som;

      /* Mark equivalent vertices and compute their contribution */

      do {

        tab_som_old_new.val[ind_som_loc] = nbr_som_new + 1;

        som_poids += 1.0;

        for (icoo = 0; icoo < 3; icoo++)
          som_coord[icoo] += (* coords)[3*ind_som_loc + icoo];

        ind_som_loc = tab_equiv_som->val[ind_som_loc];

      } while (ind_som_loc != -1);

      /* Final coordinates */

      for (icoo = 0; icoo < 3; icoo++)
        (*coords)[3*nbr_som_new + icoo] = som_coord[icoo] / som_poids;

      /* Do not forget to increment the number of vertices after merge */

      nbr_som_new += 1;
    }

  }

  /* End of main loop */
  /* ---------------- */

  /* We now know the number of vertices after merging */

  *n_vertices = nbr_som_new;
  ECS_REALLOC(*coords, nbr_som_new*3, ecs_coord_t);

  printf(_("    Number of vertices after merging  : %10d\n"), (int)nbr_som_new);

  /* Mark vertices originating from a merge */
  /*----------------------------------------*/

  nbr_som_fus = 0;

  /*
    We will not need the vertex equivalence array any more; we thus modify
    it so that only the first vertex of a linked equivalence list points
    to its first equivalent, to use it as a marker.
  */

  for (ind_som = 0; ind_som < tab_equiv_som->nbr; ind_som++) {

    if (tab_equiv_som->val[ind_som] != -1) {

      nbr_som_fus += 1;

      ind_som_loc = tab_equiv_som->val[ind_som];

      while (ind_som_loc != -1) {

        ind_som_tmp = ind_som_loc;
        ind_som_loc = tab_equiv_som->val[ind_som_loc];

        tab_equiv_som->val[ind_som_tmp] = -1;

      }
    }
  }

  printf(_("    Number of modified vertices       : %10d\n"), (int)nbr_som_fus);

  /* Return */

  return tab_som_old_new;
}

/*----------------------------------------------------------------------------
 *  Fonction qui met à jour le tableau des fusions de sommets
 *
 *  Remarque : le tableau d'équivalence (fusion) des sommets est construit de
 *             manière à ce qu'à un sommet ne comportant pas d'équivalent
 *             où de plus grand indice parmi ses équivalents corresponde la
 *             valeur -1, alors qu'un un sommet possédant des équivalents de
 *             plus grand indice corresponde le plus petit indice parmi ces
 *             équivalents (ce qui constitue une sorte de liste chaînée).
 *----------------------------------------------------------------------------*/

static void
_table_def__maj_equiv_som(size_t          ind_som_0,
                          size_t          ind_som_1,
                          ecs_tab_int_t  *tab_equiv_som)
{
  ecs_int_t  ind_som_max;
  ecs_int_t  ind_som_min;
  ecs_int_t  ind_som_tmp;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  ind_som_min = ECS_MIN(ind_som_0, ind_som_1);
  ind_som_max = ECS_MAX(ind_som_0, ind_som_1);

  /* On fusionne deux listes chaînées */
  /*----------------------------------*/

  while (   ind_som_max != ind_som_min
         && tab_equiv_som->val[ind_som_min] != ind_som_max) {

    /*
      On parcourt la liste chaînée correspondant à ind_som_min jusqu'au
      point d'insertion (si pas de point suivant dans la chaine,
      tab_equiv_som->val[ind_som_min] = -1).
    */

    while (   tab_equiv_som->val[ind_som_min] != -1
           && tab_equiv_som->val[ind_som_min] < ind_som_max)
      ind_som_min = tab_equiv_som->val[ind_som_min];

    ind_som_tmp = tab_equiv_som->val[ind_som_min];

    /*
      Si l'on est en fin de chaîne, on branche la liste chaînée correspondant
      à ind_som_max à la suite de celle correspondant à ind_som_min.
      Sinon, on doit reboucler
    */

    tab_equiv_som->val[ind_som_min] = ind_som_max;

    if (ind_som_tmp != -1) {
      ind_som_min = ind_som_max;
      ind_som_max = ind_som_tmp;
    }
  }
}

/*----------------------------------------------------------------------------
 *  Correction si nécessaire de l'orientation d'un quadrangle en connectivité
 *   nodale ; renvoie -1 si l'on ne parvient pas à corriger l'orientation
 *   ou que l'on demande une simple vérification (i.e. correc = false),
 *   0 si l'orientation initiale est bonne, et 1 en cas de correction de
 *   l'orientation par une permutation de la connectivité locale.
 *----------------------------------------------------------------------------*/

static ecs_int_t
_orient_quad(const ecs_coord_t  coord[],
             ecs_int_t          connect[],
             const bool         correc)
{
  ecs_int_t   isom_tmp;
  ecs_int_t   itri;
  ecs_int_t   nb_angle_obtus;
  ecs_coord_t  sgn;
  ecs_coord_t  normale[3];
  ecs_coord_t  vect1[3];
  ecs_coord_t  vect2[3];
  ecs_coord_t  prod_vect[3];

  ecs_int_t   passage = -1;

#define ECS_LOC_INIT_VECT(vect, i, j) ( \
  vect[0] =   coord[((connect[j-1] - 1) * 3)    ]  \
            - coord[((connect[i-1] - 1) * 3)    ], \
  vect[1] =   coord[((connect[j-1] - 1) * 3) + 1]  \
            - coord[((connect[i-1] - 1) * 3) + 1], \
  vect[2] =   coord[((connect[j-1] - 1) * 3) + 2]  \
            - coord[((connect[i-1] - 1) * 3) + 2] )

  /* Calcul d'une direction normale approchée de la face
     (peut être entrante ou sortante si la face est "croisée") */

  ECS_LOC_INIT_VECT(vect1, 1, 2);
  ECS_LOC_INIT_VECT(vect2, 1, 4);

  ECS_LOC_PRODUIT_VECTORIEL(normale, vect1, vect2);

  /* Boucle sur les renumérotations possibles */

  do {

    passage += 1;

    nb_angle_obtus = 0;

    /* Initialisation */

    switch(passage) {

    case 0:
      break;

    case 1:
      if (correc == false)
        return -1;
      isom_tmp = connect[2];
      connect[2] = connect[3];
      connect[3] = isom_tmp;
      break;

    default:
      return -1;

    }

    /* Boucle sur les coins */

    /* On compte les angles obtus, qui devraient être au nombre de 2 sur
       une face "croisée" (et 0 sur une face bien définie et convexe
       1 sur une face bien définie non convexe, soit 3 en apparaence
       sur une face bien définie convexe si la non-convexité se trouve
       au niveau du sommet 1 et que l'on a donc calculé une normale
       "à l'envers"). */

    for (itri = 0; itri < 4; itri++) {

      ECS_LOC_INIT_VECT(vect1, ((itri+2) % 4) + 1, ((itri+1) % 4) + 1);
      ECS_LOC_INIT_VECT(vect2, ( itri    % 4) + 1, ((itri+1) % 4) + 1);

      ECS_LOC_PRODUIT_VECTORIEL(prod_vect, vect1, vect2);

      /* Angle obtus si produit mixte < 0, aigu sinon. */

      sgn = ECS_LOC_PRODUIT_SCALAIRE(prod_vect, normale);

      if (sgn < 0.0)
        nb_angle_obtus += 1;

    }

  } while (nb_angle_obtus == 2);

  return (ECS_MIN(passage, 1));

#undef ECS_LOC_INIT_VECT
}

/*----------------------------------------------------------------------------
 *  Correction si nécessaire de l'orientation d'un tétraèdre en connectivité
 *   nodale ; renvoie -1 si l'on ne parvient pas à corriger l'orientation,
 *   0 si l'orientation initiale est bonne, et 1 en cas de correction de
 *   l'orientation par une permutation de la connectivité locale.
 *----------------------------------------------------------------------------*/

static ecs_int_t
_orient_tetra(const ecs_coord_t  coord[],
              ecs_int_t          connect[])
{
  ecs_int_t   isom_tmp;
  ecs_coord_t  det;
  ecs_coord_t  vect12[3];
  ecs_coord_t  vect13[3];
  ecs_coord_t  vect14[3];

  ecs_int_t   passage = -1;

#define ECS_LOC_INIT_VECT(vect, i, j) ( \
  vect[0] =   coord[((connect[j-1] - 1) * 3)    ]  \
            - coord[((connect[i-1] - 1) * 3)    ], \
  vect[1] =   coord[((connect[j-1] - 1) * 3) + 1]  \
            - coord[((connect[i-1] - 1) * 3) + 1], \
  vect[2] =   coord[((connect[j-1] - 1) * 3) + 2]  \
            - coord[((connect[i-1] - 1) * 3) + 2] )

#define ECS_LOC_PERMUTE(i, j) ( \
  isom_tmp = connect[i-1],     \
  connect[i-1] = connect[j-1], \
  connect[j-1] = isom_tmp      )

  /* Boucle sur les renumérotations possibles */

  do {

    passage += 1;

    /* Initialisation */

    switch(passage) {

    case 0:
      break;

    case 1:
      ECS_LOC_PERMUTE(2, 3);
      break;

    default: /* Retour connectivité d'origine et sortie */
      ECS_LOC_PERMUTE(2, 3);
      return -1;

    }

    ECS_LOC_INIT_VECT(vect12, 1, 2);
    ECS_LOC_INIT_VECT(vect13, 1, 3);
    ECS_LOC_INIT_VECT(vect14, 1, 4);

    det = ECS_LOC_DETERMINANT(vect12, vect13, vect14);

  } while (det < 0.0);

  return (ECS_MIN(passage, 1));

#undef ECS_LOC_INIT_VECT
#undef ECS_LOC_PERMUTE
}

/*----------------------------------------------------------------------------
 *  Correction de l'orientation issue de Foam2VTK d'un tétraèdre
 *   en connectivité nodale ; renvoie 1.
 *----------------------------------------------------------------------------*/

static ecs_int_t
_orient_tetra_of(ecs_int_t          connect[])
{
  ecs_int_t   isom_tmp;

#define ECS_LOC_PERMUTE(i, j) ( \
  isom_tmp = connect[i-1],     \
  connect[i-1] = connect[j-1], \
  connect[j-1] = isom_tmp      )

  ECS_LOC_PERMUTE(2, 3);

  return 1;

#undef ECS_LOC_PERMUTE
}

/*----------------------------------------------------------------------------
 *  Correction si nécessaire de l'orientation d'une pyramide en connectivité
 *   nodale ; renvoie -1 si l'on ne parvient pas à corriger l'orientation
 *   ou que l'on demande une simple vérification (i.e. correc = false),
 *   0 si l'orientation initiale est bonne, et 1 en cas de correction de
 *   l'orientation par une permutation de la connectivité locale.
 *
 *  Tous les cas ne sont pas détectés ou traités : on suppose que le
 *   sommet est toujours en position 5, mais que la base peut être
 *   parcourue dans l'ordre 1 4 3 2, 1 2 4 3, ou 1 3 4 2 au lieu de 1 2 3 4,
 *   dans quel cas on la réordonne.
 *----------------------------------------------------------------------------*/

static ecs_int_t
_orient_pyram(const ecs_coord_t  coord[],
              ecs_int_t          connect[],
              const bool         correc)
{
  ecs_int_t   isom;
  ecs_int_t   isom_tmp;
  ecs_int_t   connect_tmp[8];
  ecs_int_t   ret_base;
  ecs_coord_t  det;
  ecs_coord_t  vect1[3];
  ecs_coord_t  vect2[3];
  ecs_coord_t  vect3[3];

  ecs_int_t   retval = 0;
  ecs_int_t   passage = -1;

#define ECS_LOC_INIT_VECT(vect, i, j) ( \
  vect[0] =   coord[((connect_tmp[j-1] - 1) * 3)    ]  \
            - coord[((connect_tmp[i-1] - 1) * 3)    ], \
  vect[1] =   coord[((connect_tmp[j-1] - 1) * 3) + 1]  \
            - coord[((connect_tmp[i-1] - 1) * 3) + 1], \
  vect[2] =   coord[((connect_tmp[j-1] - 1) * 3) + 2]  \
            - coord[((connect_tmp[i-1] - 1) * 3) + 2] )

#define ECS_LOC_PERMUTE(i, j) ( \
  isom_tmp = connect[i-1],             \
  connect_tmp[i-1] = connect_tmp[j-1], \
  connect_tmp[j-1] = isom_tmp         )


  for (isom = 0; isom < 5; isom++)
    connect_tmp[isom] = connect[isom];

  /* Vérification et correction éventuelle de l'orientation de la base */

  ret_base = _orient_quad(coord,
                          connect_tmp,
                          correc);

  retval = ECS_MAX(ret_base, retval);

  if ((correc == false && ret_base != 0) || ret_base < 0)
    return - 1;


  /* Boucle sur les renumérotations possibles */

  do {

    passage += 1;

    /* Initialisation */

    switch(passage) {

    case 0:
      break;

    case 1:
      if (correc == false)
        return -1;
      else
        retval = 1;
      ECS_LOC_PERMUTE(2, 4);
      break;

    default: /* Retour connectivité d'origine et sortie */
      ECS_LOC_PERMUTE(2, 4);
      return -1;

    }

    ECS_LOC_INIT_VECT(vect1, 1, 2);
    ECS_LOC_INIT_VECT(vect2, 1, 4);
    ECS_LOC_INIT_VECT(vect3, 1, 5);

    det = ECS_LOC_DETERMINANT(vect1, vect2, vect3);

  } while (det < 0.0);

  if (retval > 0) {
    for (isom = 0; isom < 5; isom++)
      connect[isom] = connect_tmp[isom];
  }

  return retval;

#undef ECS_LOC_INIT_VECT
#undef ECS_LOC_PERMUTE
}

/*----------------------------------------------------------------------------
 *  Correction de l'orientation issue de Foam2VTK d'une pyramide
 *   en connectivité nodale ; renvoie 1.
 *----------------------------------------------------------------------------*/

static ecs_int_t
_orient_pyram_of(ecs_int_t          connect[])
{
  ecs_int_t   isom;
  ecs_int_t   isom_tmp;
  ecs_int_t   connect_tmp[8];

#define ECS_LOC_PERMUTE(i, j) ( \
  isom_tmp = connect[i-1],             \
  connect_tmp[i-1] = connect_tmp[j-1], \
  connect_tmp[j-1] = isom_tmp         )

  for (isom = 0; isom < 5; isom++)
    connect_tmp[isom] = connect[isom];

  /* Boucle sur les renumérotations possibles */

  ECS_LOC_PERMUTE(2, 4);

  for (isom = 0; isom < 5; isom++)
    connect[isom] = connect_tmp[isom];

  return 1;

#undef ECS_LOC_INIT_VECT
#undef ECS_LOC_PERMUTE
}

/*----------------------------------------------------------------------------
 *  Correction si nécessaire de l'orientation d'un prisme en connectivité
 *   nodale ; renvoie -1 si l'on ne parvient pas à corriger l'orientation
 *   ou que l'on demande une simple vérification (i.e. correc = false),
 *   0 si l'orientation initiale est bonne, et 1 en cas de correction de
 *   l'orientation par une permutation de la connectivité locale.
 *
 *  Tous les cas ne sont pas détectés ou traités : on suppose que les
 *   bases peuvent être parcourues dans l'ordre 1 3 2 et 4 6 5 au lieu
 *   de 1 2 3 et 4 5 6, dans quel cas on les réordonne.
 *----------------------------------------------------------------------------*/

static ecs_int_t
_orient_prism(const ecs_coord_t  coord[],
              ecs_int_t          connect[],
              const bool         correc)
{
  ecs_int_t   idim;
  ecs_int_t   isom_tmp;
  ecs_coord_t  pscal;
  ecs_coord_t  vect1[3];
  ecs_coord_t  vect2[3];
  ecs_coord_t  vect3[3];

  ecs_int_t   passage    = -1;

#define ECS_LOC_INIT_VECT(vect, i, j) ( \
  vect[0] =   coord[((connect[j-1] - 1) * 3)    ]  \
            - coord[((connect[i-1] - 1) * 3)    ], \
  vect[1] =   coord[((connect[j-1] - 1) * 3) + 1]  \
            - coord[((connect[i-1] - 1) * 3) + 1], \
  vect[2] =   coord[((connect[j-1] - 1) * 3) + 2]  \
            - coord[((connect[i-1] - 1) * 3) + 2] )

#define ECS_LOC_PERMUTE(i, j) ( \
  isom_tmp = connect[i-1],     \
  connect[i-1] = connect[j-1], \
  connect[j-1] = isom_tmp      )

  /* Boucle sur les renumérotations possibles */

  do {

    passage += 1;

    /* Initialisation */

    switch(passage) {

    case 0:
      break;

    case 1:
      if (correc == false)
        return -1;
      ECS_LOC_PERMUTE(2, 3);
      ECS_LOC_PERMUTE(5, 6);
      break;

    default: /* Retour connectivité d'origine et sortie */
      ECS_LOC_PERMUTE(2, 3);
      ECS_LOC_PERMUTE(5, 6);
      return -1;

    }

    ECS_LOC_INIT_VECT(vect2, 1, 2);
    ECS_LOC_INIT_VECT(vect3, 1, 3);

    ECS_LOC_PRODUIT_VECTORIEL(vect1, vect2, vect3);

    for (idim = 0; idim < 3; idim++)
      vect2[idim] = (  coord[((connect[4-1] - 1) * 3) + idim]
                     + coord[((connect[5-1] - 1) * 3) + idim]
                     + coord[((connect[6-1] - 1) * 3) + idim]
                     - coord[((connect[1-1] - 1) * 3) + idim]
                     - coord[((connect[2-1] - 1) * 3) + idim]
                     - coord[((connect[3-1] - 1) * 3) + idim]);

    pscal = ECS_LOC_PRODUIT_SCALAIRE(vect1, vect2);

  } while (pscal < 0.0);

  return (ECS_MIN(passage, 1));

#undef ECS_LOC_INIT_VECT
#undef ECS_LOC_PERMUTE
}

/*----------------------------------------------------------------------------
 *  Correction si nécessaire de l'orientation d'un hexaèdre en connectivité
 *   nodale ; renvoie -1 si l'on ne parvient pas à corriger l'orientation
 *   ou que l'on demande une simple vérification (i.e. correc = false),
 *   0 si l'orientation initiale est bonne, et 1 en cas de correction de
 *   l'orientation par une permutation de la connectivité locale.
 *
 *  Tous les cas ne sont pas détectés ou traités : on suppose que les
 *   bases peuvent être parcourues dans l'ordre soit 1 4 3 2 et 5 8 7 6,
 *   soit 1 2 4 3 et 5 6 8 7, soit 1 3 4 2 et 5 7 8 6 au lieu
 *   de 1 2 3 4 et 5 6 7 8, dans quel cas on les réordonne.
 *----------------------------------------------------------------------------*/

static ecs_int_t
_orient_hexa(const ecs_coord_t  coord[],
             ecs_int_t          connect[],
             const bool         correc)
{
  ecs_int_t   idim;
  ecs_int_t   isom;
  ecs_int_t   isom_tmp;
  ecs_int_t   connect_tmp[8];
  ecs_int_t   ret_base_1;
  ecs_int_t   ret_base_2;
  ecs_coord_t  pscal;
  ecs_coord_t  vect1[3];
  ecs_coord_t  vect2[3];
  ecs_coord_t  vect3[3];

  ecs_int_t   retval = 0;
  ecs_int_t   passage = -1;

#define ECS_LOC_INIT_VECT(vect, i, j) ( \
  vect[0] =   coord[((connect_tmp[j-1] - 1) * 3)    ]  \
            - coord[((connect_tmp[i-1] - 1) * 3)    ], \
  vect[1] =   coord[((connect_tmp[j-1] - 1) * 3) + 1]  \
            - coord[((connect_tmp[i-1] - 1) * 3) + 1], \
  vect[2] =   coord[((connect_tmp[j-1] - 1) * 3) + 2]  \
            - coord[((connect_tmp[i-1] - 1) * 3) + 2] )

#define ECS_LOC_PERMUTE(i, j) ( \
  isom_tmp = connect[i-1],             \
  connect_tmp[i-1] = connect_tmp[j-1], \
  connect_tmp[j-1] = isom_tmp         )


  for (isom = 0; isom < 8; isom++)
    connect_tmp[isom] = connect[isom];

  /* Vérification et correction éventuelle de l'orientation des bases */

  ret_base_1 = _orient_quad(coord,
                            connect_tmp,
                            correc);

  if ((correc == false && ret_base_1 != 0) || ret_base_1 < 0)
    return - 1;

  else if (ret_base_1 > 0) {
    ret_base_2 = _orient_quad(coord,
                              connect_tmp + 4,
                              correc);
    if (ret_base_2 != ret_base_1)
      return - 1;
    else
      retval = 1;
  }

  /* Boucle sur les renumérotations possibles */

  do {

    passage += 1;

    /* Initialisation */

    switch(passage) {

    case 0:
      break;

    case 1:
      if (correc == false)
        return -1;
      else
        retval = 1;
      ECS_LOC_PERMUTE(2, 4);
      ECS_LOC_PERMUTE(6, 8);
      break;

    default:
      return -1;
      break;

    }

    ECS_LOC_INIT_VECT(vect2, 1, 2);
    ECS_LOC_INIT_VECT(vect3, 1, 4);

    ECS_LOC_PRODUIT_VECTORIEL(vect1, vect2, vect3);

    for (idim = 0; idim < 3; idim++)
      vect2[idim] = (  coord[((connect_tmp[5-1] - 1) * 3) + idim]
                     + coord[((connect_tmp[6-1] - 1) * 3) + idim]
                     + coord[((connect_tmp[7-1] - 1) * 3) + idim]
                     + coord[((connect_tmp[8-1] - 1) * 3) + idim]
                     - coord[((connect_tmp[1-1] - 1) * 3) + idim]
                     - coord[((connect_tmp[2-1] - 1) * 3) + idim]
                     - coord[((connect_tmp[3-1] - 1) * 3) + idim]
                     - coord[((connect_tmp[4-1] - 1) * 3) + idim]) * 0.25;

    pscal = ECS_LOC_PRODUIT_SCALAIRE(vect1, vect2);

  } while (pscal < 0.0);

  if (retval > 0) {
    for (isom = 0; isom < 8; isom++)
      connect[isom] = connect_tmp[isom];
  }

  /* Vérification et correction éventuelle de l'orientation des cotés */

  connect_tmp[1-1] = connect[2-1];
  connect_tmp[2-1] = connect[3-1];
  connect_tmp[3-1] = connect[7-1];
  connect_tmp[4-1] = connect[6-1];

  ret_base_1 = _orient_quad(coord,
                            connect_tmp,
                            correc);

  if (ret_base_1 != 0)
    return - 1;

  connect_tmp[1-1] = connect[1-1];
  connect_tmp[2-1] = connect[2-1];
  connect_tmp[3-1] = connect[6-1];
  connect_tmp[4-1] = connect[5-1];

  ret_base_1 = _orient_quad(coord,
                            connect_tmp,
                            correc);

  if (ret_base_1 != 0)
    return - 1;

  return retval;

#undef ECS_LOC_INIT_VECT
#undef ECS_LOC_PERMUTE
}

/*----------------------------------------------------------------------------
 *  Correction de l'orientation issue de Foam2VTK d'un hexaedre
 *   en connectivité nodale ; renvoie 1.
 *----------------------------------------------------------------------------*/

static ecs_int_t
_orient_hexa_of(ecs_int_t          connect[])
{
  ecs_int_t   isom;
  ecs_int_t   isom_tmp;
  ecs_int_t   connect_tmp[8];

#define ECS_LOC_PERMUTE(i, j) ( \
  isom_tmp = connect[i-1],             \
  connect_tmp[i-1] = connect_tmp[j-1], \
  connect_tmp[j-1] = isom_tmp         )

  for (isom = 0; isom < 8; isom++)
    connect_tmp[isom] = connect[isom];

  /* Renumbering */

  ECS_LOC_PERMUTE(2, 4);
  ECS_LOC_PERMUTE(6, 8);

  for (isom = 0; isom < 8; isom++)
    connect[isom] = connect_tmp[isom];

  return 1;

#undef ECS_LOC_PERMUTE
}

/*----------------------------------------------------------------------------
 * Compute contribution of a polygonal face to a cell's volume,
 *  using Stoke's theorem
 *
 *                          Pi+1
 *              *---------*                   B  : barycentre of the polygon
 *             / .       . \
 *            /   .     .   \                 Pi : vertices of the polygon
 *           /     .   .     \
 *          /       . .  Ti   \               Ti : triangle
 *         *.........B.........* Pi
 *     Pn-1 \       . .       /
 *           \     .   .     /
 *            \   .     .   /
 *             \ .   T0  . /
 *              *---------*
 *            P0
 *----------------------------------------------------------------------------*/

static double
_compute_face_vol_contrib(ecs_int_t          n_face_vertices,
                          const ecs_coord_t  vtx_coord[],
                          const ecs_int_t    face_vtx_lst[])
{
  ecs_int_t   i, tri_id, vtx_id;
  ecs_coord_t  tri_center[3], tri_norm[3];
  ecs_coord_t  vect1[3], vect2[3];

  ecs_coord_t  face_barycentre[3] = {0., 0., 0.};

  double inv3 = 1./3.;

  double vol_contrib = 0.;

  /* Compute barycentre (B) coordinates for the polygon (P) */

  for (vtx_id = 0; vtx_id < n_face_vertices; vtx_id++) {
    size_t coord_id = face_vtx_lst[vtx_id] - 1;
    for (i = 0; i < 3; i++)
      face_barycentre[i] += vtx_coord[coord_id*3 + i];
  }

  for (i = 0; i < 3; i++)
    face_barycentre[i] /= n_face_vertices;

  /* Loop on triangles of the face */

  for (tri_id = 0; tri_id < n_face_vertices; tri_id++) {

    size_t coord_id_1 = face_vtx_lst[tri_id % n_face_vertices] - 1;
    size_t coord_id_2 = face_vtx_lst[(tri_id + 1)% n_face_vertices] - 1;

    for (i = 0; i < 3; i++) {
      vect1[i] = vtx_coord[coord_id_1*3 + i] - face_barycentre[i];
      vect2[i] = vtx_coord[coord_id_2*3 + i] - face_barycentre[i];
    }

    ECS_LOC_PRODUIT_VECTORIEL(tri_norm, vect1, vect2);

    for (i = 0; i < 3; i++)
      tri_norm[i] *= 0.5;

    /* Computation of the center of the triangle Ti */

    for (i = 0; i < 3; i++)
      tri_center[i] = (  face_barycentre[i]
                       + vtx_coord[coord_id_1*3 + i]
                       + vtx_coord[coord_id_2*3 + i]) * inv3;

    /* Contribution to cell volume using Stoke's formula */

    for (i = 0; i < 3; i++)
      vol_contrib += (tri_norm[i] * tri_center[i]);

  } /* End of loop on triangles of the face */

  return vol_contrib;
}

/*----------------------------------------------------------------------------
 *  Correction si nécessaire de l'orientation d'un polyedre en connectivité
 *   nodale ; renvoie -1 si l'on ne parvient pas à corriger l'orientation
 *   ou que l'on demande une simple vérification (i.e. correc = false),
 *   0 si l'orientation initiale est bonne, et 1 en cas de correction de
 *   l'orientation par une permutation de la connectivité locale.
 *
 *  On suppose que toutes les faces polygonales sont bien définies, mais
 *   que certaines faces peuvent être définies avec une normale intérieure
 *   à la cellule et non pas extérieure. On suppose que la non-convexité
 *   du polyèdre est limitée.
 *----------------------------------------------------------------------------*/

static ecs_int_t
_orient_polyhedron(const ecs_coord_t   coord[],
                   ecs_int_t           connect[],
                   ecs_int_t           size,
                   bool                correc,
                   ecs_tab_int_t      *face_index,
                   ecs_tab_int_t      *face_marker,
                   ecs_tab_int_t      *edges)
{
  size_t      i, j, face_id, edge_id;

  double      cell_vol = 0.;
  size_t      n_unmarked_faces = 0;
  size_t      n_faces = 0;
  size_t      n_edges = 0;
  ecs_int_t   retval = 0;

  /* Ensure working arrays are large enough */

  /* Mark faces */

  {
    ecs_int_t   premier_som  = -1;

    for (i = 0; i < (size_t)size; i++) {

      if (face_index->nbr < n_faces + 1) {
        face_index->nbr = ECS_MAX(n_faces + 1, face_index->nbr*2);
        ECS_REALLOC(face_index->val, face_index->nbr, ecs_int_t);
      }

      if (premier_som == -1) {
        face_index->val[n_faces] = i;
        premier_som = connect[i];
      }
      else if (connect[i] == premier_som) {
        n_faces += 1;
        premier_som = -1;
      }
    }

    if (face_index->nbr < n_faces + 1) {
      face_index->nbr = ECS_MAX(n_faces + 1, face_index->nbr*2);
      ECS_REALLOC(face_index->val, face_index->nbr, ecs_int_t);
    }

    face_index->val[n_faces] = size;

    /* face_marker: 0 initially, 1 if oriented, -1 if inverted */

    if (face_marker->nbr < n_faces) {
      face_marker->nbr = ECS_MAX(n_faces, face_marker->nbr*2);
      ECS_REALLOC(face_marker->val, face_marker->nbr, ecs_int_t);
    }

    for (face_id = 0; face_id < n_faces; face_id++)
      face_marker->val[face_id] = 0;
  }

  if (n_faces == 0)
    return 0;

  /* Extract edges and build edge-face connectivity */
  /*------------------------------------------------*/

  for (face_id = 0; face_id < n_faces; face_id++) {

    size_t start_id = face_index->val[face_id];
    size_t end_id = face_index->val[face_id + 1] - 1;

    /* Loop on unassociated edges */

    for (i = start_id; i < end_id; i++) {

      ecs_int_t v_id_1 = connect[i], v_id_2 = connect[i + 1];
      ecs_int_t is_match = 0;

      /* Compare with other edges */

      for (edge_id = 0; edge_id < n_edges; edge_id++) {

        ecs_int_t v_id_1_cmp = edges->val[edge_id*4];
        ecs_int_t v_id_2_cmp = edges->val[edge_id*4 + 1];

        if ((v_id_1 == v_id_2_cmp) && (v_id_2 == v_id_1_cmp))
          is_match = 1;

        else if ((v_id_1 == v_id_1_cmp) && (v_id_2 == v_id_2_cmp))
          is_match = -1;

        /* Each edge should be shared by exactly 2 faces */

        if (is_match != 0) {
          if (edges->val[edge_id*4 + 3] == 0) {
            edges->val[edge_id*4 + 3] = is_match * (face_id + 1);
            break;
          }
          else
            return -1;
        }
      }

      /* If edge is unmatched, add it */

      if (is_match == 0) {

        if (edges->nbr < (n_edges + 1)*4) {
          edges->nbr = ECS_MAX(edges->nbr*2, (n_edges + 1)*4);
          ECS_REALLOC(edges->val, edges->nbr, ecs_int_t);
        }

        edges->val[n_edges*4] = v_id_1;
        edges->val[n_edges*4 + 1] = v_id_2;
        edges->val[n_edges*4 + 2] = face_id + 1;
        edges->val[n_edges*4 + 3] = 0;

        n_edges += 1;
      }

    }

  }

  /* Check if each edge is associated with 2 faces
     (i.e., that cell is closed) */

  for (edge_id = 0; edge_id < n_edges; edge_id++) {
    if (edges->val[edge_id*4 + 3] == 0)
      return -1;
  }

  /* Now check if face orientations are compatible */
  /*-----------------------------------------------*/

  face_marker->val[0] = 1; /* First face defines reference */
  n_unmarked_faces = n_faces - 1;

  while (n_unmarked_faces > 0) {

    for (edge_id = 0; edge_id < n_edges; edge_id++) {

      ecs_int_t face_num_2 = edges->val[edge_id*4 + 3];
      ecs_int_t face_id_1 = edges->val[edge_id*4 + 2] - 1;
      ecs_int_t face_id_2 = ECS_ABS(face_num_2) - 1;

      if (   face_marker->val[face_id_1] == 0
          && face_marker->val[face_id_2] != 0) {
        if (face_num_2 > 0)
          face_marker->val[face_id_1] = face_marker->val[face_id_2];
        else
          face_marker->val[face_id_1] = - face_marker->val[face_id_2];
        n_unmarked_faces -= 1;
      }

      else if (   face_marker->val[face_id_1] != 0
               && face_marker->val[face_id_2] == 0) {
        if (face_num_2 > 0)
          face_marker->val[face_id_2] = face_marker->val[face_id_1];
        else
          face_marker->val[face_id_2] = - face_marker->val[face_id_1];
        n_unmarked_faces -= 1;
      }
    }

  }

  for (edge_id = 0; edge_id < n_edges; edge_id++) {

    ecs_int_t face_num_2 = edges->val[edge_id*4 + 3];
    ecs_int_t face_id_1 = edges->val[edge_id*4 + 2] - 1;
    ecs_int_t face_id_2 = ECS_ABS(face_num_2) - 1;

    ecs_int_t marker_product
      = face_marker->val[face_id_1]*face_marker->val[face_id_2];

    if (   (face_num_2 < 0 && marker_product > 0)
        || (face_num_2 > 0 && marker_product < 0))
      return -1;

  }

  /* At this stage, topology is correct; check for inside-out cell */
  /*---------------------------------------------------------------*/

  for (face_id = 0; face_id < n_faces; face_id++) {

    size_t start_id = face_index->val[face_id];
    size_t end_id = face_index->val[face_id + 1] - 1;

    double vol_contrib
      = _compute_face_vol_contrib(end_id - start_id,
                                  coord,
                                  connect + start_id);
    cell_vol += vol_contrib * (double)(face_marker->val[face_id]);
  }

  /* Invert orientation if cell is inside_out */

  if (cell_vol < 0.) {
    for (face_id = 0; face_id < n_faces; face_id++)
      face_marker->val[face_id] *= -1;
  }

  /* Now correct connectivity if required */
  /*--------------------------------------*/

  if (correc == true) {

    for (face_id = 0; face_id < n_faces; face_id++) {

      if (face_marker->val[face_id] == -1) {

        for (i = face_index->val[face_id] + 1,
               j = face_index->val[face_id + 1] - 2;
             i < j;
             i++, j--) {
          ecs_int_t connect_tmp = connect[i];
          connect[i] = connect[j];
          connect[j] = connect_tmp;
        }

        retval = 1;

      }
    }
  }

  else {

    for (face_id = 0; face_id < n_faces; face_id++) {
      if (face_marker->val[face_id] == -1)
        retval = -1;
    }
  }

  return retval;
}

/*============================================================================
 *                             Fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction qui réalise le tri des types géométriques
 *  La fonction affiche le nombre d'éléments par type géométrique
 *----------------------------------------------------------------------------*/

ecs_tab_int_t
ecs_table_def__trie_typ(ecs_table_t  *this_table_def,
                        int           dim_elt)
{
  size_t       ielt;
  size_t       nbr_elt;
  size_t       nbr_som;

  ecs_tab_int_t   tab_typ_geo_ord;
  ecs_tab_int_t   tab_typ_geo;
  ecs_tab_int_t   vect_renum;

  ecs_elt_typ_t   typ_geo;

  const ecs_elt_typ_t  *typ_geo_base;

  const int lng_imp = ECS_LNG_AFF_STR - ECS_LOC_LNG_MAX_NOM_TYP;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(this_table_def != NULL);

  vect_renum.nbr = 0   ;
  vect_renum.val = NULL;

  nbr_elt = this_table_def->nbr;

  if (nbr_elt == 0)
    return vect_renum;

  /* Détermination de base du type d'élément "classique" */

  assert(dim_elt >=2 && dim_elt <= 3);

  typ_geo_base = ecs_glob_typ_elt[dim_elt - 2];

  /* Si tous les éléments sont de même type, rien à faire */
  /*------------------------------------------------------*/

  if (this_table_def->pos == NULL) {

    if (this_table_def->pas < 9)
      typ_geo = typ_geo_base[this_table_def->pas];
    else if (dim_elt == 2)
      typ_geo = ECS_ELT_TYP_FAC_POLY;
    else if (dim_elt == 3)
      typ_geo = ECS_ELT_TYP_CEL_POLY;
    else
      typ_geo = ECS_ELT_TYP_NUL;

    printf("  %-*s%*.*s : %*lu\n",
           lng_imp, _("Number of elements"),
           ECS_LOC_LNG_MAX_NOM_TYP, ECS_LOC_LNG_MAX_NOM_TYP,
           ecs_fic_elt_typ_liste_c[typ_geo].nom,
           ECS_LNG_AFF_ENT, (unsigned long)nbr_elt);

  }

  /* Si tous les éléments ne sont pas de même type */
  /*-----------------------------------------------*/

  else {

    /* Construction d'un tableau temporaire de type géométrique */

    tab_typ_geo.nbr = nbr_elt;
    ECS_MALLOC(tab_typ_geo.val, tab_typ_geo.nbr, ecs_int_t);

    for (ielt = 0; ielt < nbr_elt; ielt++) {

      nbr_som =   this_table_def->pos[ielt + 1]
                - this_table_def->pos[ielt];

      if (nbr_som < 9)
        tab_typ_geo.val[ielt] = typ_geo_base[nbr_som];

      else if (dim_elt == 2)
        tab_typ_geo.val[ielt] = ECS_ELT_TYP_FAC_POLY;

      else if (dim_elt == 3)
        tab_typ_geo.val[ielt] = ECS_ELT_TYP_CEL_POLY;

      else
        tab_typ_geo.val[ielt] = ECS_ELT_TYP_NUL;

    }

    /* On regarde si les types géométriques ne sont pas deja ordonnés */

    ielt = 1;
    while (ielt < nbr_elt                                     &&
           tab_typ_geo.val[ielt] >= tab_typ_geo.val[ielt - 1]   )
      ielt++;

    if (ielt < nbr_elt) {

      /* Les types géométriques ne sont pas ordonnés */
      /* On ordonne les types géométriques */

      vect_renum.nbr = nbr_elt;
      ECS_MALLOC(vect_renum.val, nbr_elt, ecs_int_t);


      tab_typ_geo_ord = ecs_tab_int__trie_et_renvoie(tab_typ_geo,
                                                     vect_renum);

      /*
        `vect_renum' prend pour indice les indices nouveaux, et ses valeurs
        contiennent les indices anciens correspondants
        On inverse le contenu de `vect_renum' :
        à chaque indice ancien, `vect_renum' donne la valeur du nouvel indice
      */

      ecs_tab_int__inverse(&vect_renum);

      ECS_FREE(tab_typ_geo.val);

      tab_typ_geo.val = tab_typ_geo_ord.val;

      tab_typ_geo_ord.nbr = 0;
      tab_typ_geo_ord.val = NULL;

    }

    /* Message d'information sur la composition des éléments du maillage */

    {
      ecs_int_t val_typ_ref;
      size_t    nbr_val_typ_geo;

      size_t    cpt_ielt = 0;

      while (cpt_ielt < nbr_elt) {

        val_typ_ref = tab_typ_geo.val[cpt_ielt];

        /* On compte le nombre d'éléments ayant le même type géométrique */

        nbr_val_typ_geo = 0;

        for (ielt = cpt_ielt;
             ielt < nbr_elt && tab_typ_geo.val[ielt] == val_typ_ref; ielt++)
          nbr_val_typ_geo++;

        printf("  %-*s%*.*s : %*lu\n",
               lng_imp, _("Number of elements"),
               ECS_LOC_LNG_MAX_NOM_TYP, ECS_LOC_LNG_MAX_NOM_TYP,
               ecs_fic_elt_typ_liste_c[val_typ_ref].nom,
               ECS_LNG_AFF_ENT, (unsigned long)nbr_val_typ_geo);

        cpt_ielt += nbr_val_typ_geo;

      }
    }

    /* Le tableau de type géométrique n'est plus nécessaire */

    tab_typ_geo.nbr = 0;
    ECS_FREE(tab_typ_geo.val);

  }

  /* Renvoi du vecteur de renumérotation */
  /*-------------------------------------*/

  return vect_renum;
}

/*----------------------------------------------------------------------------
 *  Fonction qui construit
 *   les définitions des faces par décomposition des tables des cellules
 *----------------------------------------------------------------------------*/

void
ecs_table_def__decompose_cel(ecs_table_t  **table_def_fac,
                             ecs_table_t   *table_def_cel)
{
  size_t      nbr_cel;
  size_t      nbr_def;
  size_t      nbr_fac;
  size_t      nbr_fac_old;
  size_t      nbr_val_cel;
  size_t      nbr_val_fac;
  size_t      nbr_val_fac_old;
  size_t      ind_pos_cel;
  size_t      ind_pos_loc;
  size_t      ind_pos_sui;
  size_t      nbr_pos_loc;
  ecs_int_t   num_def;
  ecs_int_t   marqueur_fin;

  size_t      isom;
  size_t      icel;
  int         idef;
  size_t      cpt_fac;
  int         ifac;

  ecs_int_t   typ_geo_cel;        /* Type géométrique élément en cours */

  ecs_int_t typ_geo_base[9] = {ECS_ELT_TYP_NUL,
                               ECS_ELT_TYP_NUL,
                               ECS_ELT_TYP_NUL,
                               ECS_ELT_TYP_NUL,
                               ECS_ELT_TYP_CEL_TETRA,
                               ECS_ELT_TYP_CEL_PYRAM,
                               ECS_ELT_TYP_CEL_PRISM,
                               ECS_ELT_TYP_NUL,
                               ECS_ELT_TYP_CEL_HEXA};

  ecs_size_t   *def_cel_fac_pos = NULL;
  ecs_int_t    *def_cel_fac_val = NULL;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /*=================*/
  /* Initialisations */
  /*=================*/

  nbr_cel  = table_def_cel->nbr;

  ecs_table__regle_en_pos(table_def_cel);

  /* Construction, pour les cellules, des tables principaux */
  /*--------------------------------------------------------*/

  /* Boucle de comptage pour l'allocation des sous-éléments */
  /*--------------------------------------------------------*/

  if (*table_def_fac != NULL) {
    nbr_fac_old = (*table_def_fac)->nbr;
    nbr_val_fac_old = ecs_table__ret_val_nbr(*table_def_fac);
  }
  else {
    nbr_fac_old = 0;
    nbr_val_fac_old = 0;
  }

  nbr_fac     = 0;
  nbr_val_fac = 0;
  nbr_def     = 0;

  for (icel = 0; icel < nbr_cel; icel++ ) {

    ind_pos_cel = table_def_cel->pos[icel] - 1;

    nbr_val_cel = table_def_cel->pos[icel + 1] - 1 - ind_pos_cel;

    /* Traitement des cellules "classiques" */
    /*--------------------------------------*/

    if (nbr_val_cel < 9) {

      typ_geo_cel = typ_geo_base[nbr_val_cel];

      /* Boucle sur les sous-éléments définissant la cellulle */

      for (ifac = 0;
           ifac < ecs_fic_elt_typ_liste_c[typ_geo_cel].nbr_sous_elt;
           ifac++ ) {

        const ecs_sous_elt_t  * sous_elt
          = &(ecs_fic_elt_typ_liste_c[typ_geo_cel].sous_elt[ifac]);

        nbr_fac++;

        for (idef = 0;
             idef < (ecs_fic_elt_typ_liste_c[sous_elt->elt_typ]).nbr_som;
             idef++);

        nbr_val_fac += idef;

      }

    }

    /* Traitement des éléments de type "polyèdre" */
    /*--------------------------------------------*/

    else {

      ind_pos_sui = table_def_cel->pos[icel + 1] - 1;
      nbr_pos_loc = ind_pos_sui - ind_pos_cel;

      /* Convention : définition nodale cellule->sommets avec numéros de
         premiers sommets répétés en fin de liste pour marquer la fin
         de chaque face */

      marqueur_fin = -1;

      for (isom = ind_pos_cel; isom < ind_pos_sui; isom++) {

        if (table_def_cel->val[isom] != marqueur_fin) {
          nbr_val_fac += 1;
          if (marqueur_fin == -1)
            marqueur_fin = table_def_cel->val[isom];
        }
        else {
          marqueur_fin = -1;
          nbr_fac += 1;
        }

      }

    }

  } /* Fin de la boucle de comptage sur les cellules */

  /* Allocation et initialisation pour les faces  */
  /*  des tableaux associés aux définitions       */
  /*----------------------------------------------*/

  nbr_fac += nbr_fac_old;
  nbr_val_fac += nbr_val_fac_old;

  if (*table_def_fac != NULL) {
    ecs_table__regle_en_pos(*table_def_fac);
    (*table_def_fac)->nbr = nbr_fac;
    ECS_REALLOC((*table_def_fac)->pos, nbr_fac + 1, ecs_size_t);
    ECS_REALLOC((*table_def_fac)->val, nbr_val_fac, ecs_int_t);
  }
  else {
    *table_def_fac = ecs_table__alloue(nbr_fac,
                                       nbr_val_fac);
    (*table_def_fac)->pos[0] = 1;
  }

  ECS_MALLOC(def_cel_fac_pos, nbr_cel + 1, ecs_size_t);
  ECS_MALLOC(def_cel_fac_val, nbr_fac - nbr_fac_old, ecs_int_t);

  def_cel_fac_pos[0] = 1;

  /*=======================================*/
  /* Boucle sur les cellules à transformer */
  /*=======================================*/

  cpt_fac = nbr_fac_old;
  nbr_def = nbr_val_fac_old;

  for (icel = 0; icel < nbr_cel; icel++ ) {

    ind_pos_cel = table_def_cel->pos[icel] - 1;

    nbr_val_cel = table_def_cel->pos[icel + 1] - 1 - ind_pos_cel;

    /*--------------------------------------*/
    /* Traitement des éléments "classiques" */
    /*--------------------------------------*/

    if (nbr_val_cel < 9) {

      typ_geo_cel = typ_geo_base[nbr_val_cel];

      /* Boucle sur les faces définissant la cellulle */
      /*==============================================*/

      for (ifac = 0;
           ifac < ecs_fic_elt_typ_liste_c[typ_geo_cel].nbr_sous_elt;
           ifac++ ) {

        const ecs_sous_elt_t  * sous_elt
          = &(ecs_fic_elt_typ_liste_c[typ_geo_cel].sous_elt[ifac]);

        /* Boucle sur les sommets définissant la face */
        /*--------------------------------------------*/

        for (idef = 0;
             idef < (ecs_fic_elt_typ_liste_c[sous_elt->elt_typ]).nbr_som;
             idef++) {

          /* Définition de la face en fonction des sommets */

          num_def = sous_elt->som[idef];

          (*table_def_fac)->val[nbr_def++]
            = table_def_cel->val[ind_pos_cel + num_def - 1];

        }

        /* Position de la face dans sa définition en fonction des sommets */
        /*----------------------------------------------------------------*/

        (*table_def_fac)->pos[cpt_fac + 1]
          = (*table_def_fac)->pos[cpt_fac] + idef;

        /* Détermination de la cellule en fonction des faces */
        /*---------------------------------------------------*/

        def_cel_fac_val[cpt_fac - nbr_fac_old] = cpt_fac + 1;

        cpt_fac++;

      } /* Fin de la boucle sur les faces d'une cellule */

    }


    /*--------------------------------------------*/
    /* Traitement des éléments de type "polyèdre" */
    /*--------------------------------------------*/

    else {

      /* Boucle sur les faces définissant le polyèdre */
      /*==============================================*/

      ifac = 0;
      marqueur_fin = -1;

      nbr_pos_loc = table_def_cel->pos[icel + 1] - 1 - ind_pos_cel;

      for (ind_pos_loc = 0; ind_pos_loc < nbr_pos_loc; ind_pos_loc++) {

        /* Définition de la face en fonction des sommets */

        if (table_def_cel->val[ind_pos_cel + ind_pos_loc] != marqueur_fin) {

          (*table_def_fac)->val[nbr_def++]
            = table_def_cel->val[ind_pos_cel + ind_pos_loc];

          if (marqueur_fin == -1)
            marqueur_fin = table_def_cel->val[ind_pos_cel + ind_pos_loc];

        }

        /* Position de la face dans sa définition en fonction des sommets */

        else {

          (*table_def_fac)->pos[cpt_fac + 1] = nbr_def + 1;

          marqueur_fin = -1;

          def_cel_fac_val[cpt_fac - nbr_fac_old] = cpt_fac + 1;

          /* Incrémentation du nombre de faces */

          ifac++;
          cpt_fac++;

        }

      } /* Fin de la boucle sur les faces du polyèdre */

    }

    /* A ce point, on a : ifac
       = ecs_fic_elt_typ_liste_c[table_typ_geo_cel->val[icel]].nbr_sous_elt
       pour des éléments classiques, et ifac est égal au nombre de faces pour
       des éléments polyédriques. */

    def_cel_fac_pos[icel + 1] = def_cel_fac_pos[icel] + ifac;

  } /* Fin de la boucle sur les éléments */

  /* Mise à jour de la table */
  /*-------------------------*/

  ECS_FREE(table_def_cel->pos);
  ECS_FREE(table_def_cel->val);

  table_def_cel->pos = def_cel_fac_pos;
  table_def_cel->val = def_cel_fac_val;

  /* ecs_table__pos_en_regle(*table_def_fac); */
  ecs_table__pos_en_regle(table_def_cel);
}

/*----------------------------------------------------------------------------
 *  Fonction qui realise la fusion des definitions des elements
 *----------------------------------------------------------------------------*/

ecs_tab_int_t
ecs_table_def__fusionne(ecs_table_t    *table_def,
                        size_t         *nbr_elt_cpct,
                        ecs_tab_int_t  *signe_elt)
{
  size_t           cpt_sup_fin;
  size_t           ind_cmp;
  size_t           ind_inf;
  size_t           ind_loc_sup;
  size_t           ind_loc_cmp;
  size_t           ind_pos;
  size_t           ind_pos_loc;
  size_t           ind_sup;
  size_t           nbr_sup_ini;
  size_t           nbr_inf;
  ecs_int_t        num_inf;
  ecs_int_t        num_inf_loc;
  ecs_int_t        num_inf_min;
  ecs_int_t        num_inf_min_cmp;
  size_t           pos_cmp;
  size_t           pos_cpt;
  size_t           pos_sup;
  int              sgn;

  size_t           ind_pos_sup[3];
  size_t           ind_pos_cmp[3];

  ecs_tab_int_t    cpt_ref_inf;

  ecs_size_t      *pos_recherche = NULL;
  ecs_int_t       *val_recherche = NULL;

  ecs_tab_int_t    tab_transf;    /* Tableau de transformation */

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  ecs_table__regle_en_pos(table_def);

  /* Fusion des definitions des elements sur la liste initiale  */
  /*------------------------------------------------------------*/

  /* Initialisations avec précautions pour cas vide */

  tab_transf.nbr = 0;
  tab_transf.val = NULL;

  if (table_def == NULL)
    return tab_transf;

  if (table_def == NULL)
    return tab_transf;

  nbr_sup_ini = table_def->nbr;
  cpt_sup_fin = 0;

  if (nbr_sup_ini < 1)
    return tab_transf;

  if (signe_elt != NULL) {
    signe_elt->nbr = nbr_sup_ini;
    ECS_MALLOC(signe_elt->val, nbr_sup_ini, ecs_int_t);
  }

  /* Comptage du nombre de sous-entités */

  nbr_inf = 0;

  for (ind_sup = 0;
       ind_sup < table_def->pos[nbr_sup_ini] - 1;
       ind_sup++) {
    num_inf = table_def->val[ind_sup];
    if ((size_t)(ECS_ABS(num_inf)) > nbr_inf)
      nbr_inf = ECS_ABS(num_inf);
  }

  /* Construction du tableau de recherche des sous-entités */
  /*-------------------------------------------------------*/

  ECS_MALLOC(pos_recherche, nbr_inf + 1, ecs_size_t);
  ECS_MALLOC(val_recherche, table_def->nbr, ecs_int_t);

  /*
    Comptage pour chaque élément de l'entité inférieure du nombre de fois où il
    apparaît comme élement d'indice minimal d'un élément de l'entité supérieure
  */

  cpt_ref_inf = ecs_tab_int__cree_init(nbr_inf, 0);

  for (ind_sup = 0; ind_sup < table_def->nbr; ind_sup++) {

    ind_pos_loc = table_def->pos[ind_sup] - 1;

    num_inf_min = ECS_ABS(table_def->val[ind_pos_loc]);
    while (++ind_pos_loc < table_def->pos[ind_sup + 1] -1) {
      num_inf_loc = ECS_ABS(table_def->val[ind_pos_loc]);
      if (num_inf_loc < num_inf_min)
        num_inf_min = num_inf_loc;
    }

    assert(num_inf_min > 0 && (size_t)num_inf_min <= nbr_inf);

    cpt_ref_inf.val[num_inf_min - 1] += 1;
  }


  /* Construction du vecteur */

  pos_recherche[0] = 1;

  for (ind_inf = 0; ind_inf < nbr_inf; ind_inf++) {

    pos_recherche[ind_inf + 1]
      = pos_recherche[ind_inf] + cpt_ref_inf.val[ind_inf];

    cpt_ref_inf.val[ind_inf] = 0;

  }


  for (ind_sup = 0; ind_sup < table_def->nbr; ind_sup++) {

    ind_pos_loc = table_def->pos[ind_sup] - 1;

    num_inf_min = ECS_ABS(table_def->val[ind_pos_loc]);
    while (ind_pos_loc < table_def->pos[ind_sup + 1] -1) {
      num_inf_loc = ECS_ABS(table_def->val[ind_pos_loc]);
      ind_pos_loc++;
      if (num_inf_loc < num_inf_min)
        num_inf_min = num_inf_loc;
    }

    ind_inf = num_inf_min - 1;

    ind_pos =   pos_recherche[ind_inf] - 1
              + cpt_ref_inf.val[ind_inf];

    cpt_ref_inf.val[ind_inf] += 1;

    val_recherche[ind_pos] = ind_sup + 1;

  }


  /* Libération du tableau auxiliaire */

  cpt_ref_inf.nbr = 0;
  ECS_FREE(cpt_ref_inf.val);

  /* Allocation du tableau de transformation */

  tab_transf = ecs_tab_int__cree(nbr_sup_ini);

  for (ind_sup = 0; ind_sup < nbr_sup_ini; ind_sup++)
    tab_transf.val[ind_sup] = ind_sup;


  /* Boucle principale de recherche sur les éléments supérieurs */
  /*------------------------------------------------------------*/

  /*
    Le premier élément ne peut pas être fusionné avec un élément précédent,
    la boucle commence donc au deuxième
  */

  tab_transf.val[0] = 0;

  if (signe_elt != NULL)
    signe_elt->val[0] = 1;

  cpt_sup_fin       = 1;

  for (ind_sup = 1; ind_sup < table_def->nbr; ind_sup++) {

    /* Recherche élement entité inférieure de plus petit numéro référencé */

    ind_pos_sup[0] = table_def->pos[ind_sup    ] - 1; /* début */
    ind_pos_sup[1] = table_def->pos[ind_sup + 1] - 1; /* fin */
    ind_pos_sup[2] = table_def->pos[ind_sup    ] - 1; /* plus petit */

    ind_pos_loc = ind_pos_sup[0];

    num_inf_min = ECS_ABS(table_def->val[ind_pos_loc]);
    while (++ind_pos_loc < ind_pos_sup[1]) {
      num_inf_loc = ECS_ABS(table_def->val[ind_pos_loc]);
      if (num_inf_loc < num_inf_min) {
        num_inf_min    = num_inf_loc;
        ind_pos_sup[2] = ind_pos_loc;
      }
    }

    /*
      On cherche des éléments de l'entité courante de plus petit numéro que
      l'entité courante ayant même plus petit élément de l'entité inférieure
      (recherche de candidats pour la fusion)
    */

    ind_inf = num_inf_min - 1;
    sgn     = 0;

    for (pos_cmp = pos_recherche[ind_inf]     - 1;
         pos_cmp < pos_recherche[ind_inf + 1] - 1;
         pos_cmp++) {

      ind_cmp = val_recherche[pos_cmp] - 1;

      /* Repérage point de départ pour comparaison */

      if (ind_cmp < ind_sup) {

        ind_pos_cmp[0] = table_def->pos[ind_cmp    ] - 1; /* début */
        ind_pos_cmp[1] = table_def->pos[ind_cmp + 1] - 1; /* fin */
        ind_pos_cmp[2] = table_def->pos[ind_cmp    ] - 1; /* plus petit */

        assert(ind_pos_cmp[1] > ind_pos_cmp[0]);

        ind_pos_cmp[1] = table_def->pos[ind_cmp + 1] - 1;  /* fin */

        ind_pos_loc = ind_pos_cmp[0];

        num_inf_min_cmp = ECS_ABS(table_def->val[ind_pos_loc]);
        while ((++ind_pos_loc) < ind_pos_cmp[1]) {
          num_inf_loc = ECS_ABS(table_def->val[ind_pos_loc]);
          if (num_inf_loc < num_inf_min_cmp) {
            num_inf_min_cmp = num_inf_loc;
            ind_pos_cmp[2]  = ind_pos_loc;
          }
        }

        /* Comparaison des définitions */

        for (sgn = 1; sgn > -2; sgn -= 2) {

          ind_loc_sup = ind_pos_sup[2];
          ind_loc_cmp = ind_pos_cmp[2];

          do {

            ind_loc_sup++;
            if (ind_loc_sup == ind_pos_sup[1])
              ind_loc_sup = ind_pos_sup[0];

            ind_loc_cmp += sgn;
            if (ind_loc_cmp == ind_pos_cmp[1])
              ind_loc_cmp = ind_pos_cmp[0];
            else if (   ind_loc_cmp < ind_pos_cmp[0]
                     || ind_loc_cmp > ind_pos_cmp[1])
              ind_loc_cmp = ind_pos_cmp[1] - 1;

          } while (   (   ECS_ABS(table_def->val[ind_loc_sup])
                       == ECS_ABS(table_def->val[ind_loc_cmp]))
                   && ind_loc_sup != ind_pos_sup[2]
                   && ind_loc_cmp != ind_pos_cmp[2]);

          if (   ind_loc_sup == ind_pos_sup[2]
              && ind_loc_cmp == ind_pos_cmp[2])
            break; /* Sortie boucle sur signe parcours (1, -1, -3) */

        }

        /*
          Si sgn =  1, les entités sont confondues, de même sens;
          Si sgn = -1, elles sont confondues, de sens inverse;
          Sinon, elles ne sont pas confondues
        */

        if (sgn == 1 || sgn == -1)
          break; /* Sortie boucle sur pos_cmp pour recherche de candidats */

      } /* Fin de la comparaison référence/candidat (if ind_cmp < ind_sup) */

    } /* Fin boucle sur pos_cmp recherche de candidats */


    /* Si on a trouvé une entité à fusionner */

    if (sgn == 1 || sgn == -1) {

      assert(pos_cmp < pos_recherche[ind_inf + 1] - 1);

      if (sgn == 1 && (ind_pos_cmp[1] - ind_pos_cmp[0] == 2)) {

        /*
          Si l'on a deux entités inférieures par entité supérieure,
          la permutation cyclique peut trouver une égalité quel
          que soit le signe; on le corrige si nécessaire
        */

        if (   ECS_ABS(table_def->val[ind_pos_sup[0]])
            != ECS_ABS(table_def->val[ind_pos_cmp[0]]))
          sgn = -1;

      }

      tab_transf.val[ind_sup] = tab_transf.val[ind_cmp];

      if (signe_elt != NULL)
        signe_elt->val[ind_sup] = sgn;

    }
    else {

      tab_transf.val[ind_sup] = cpt_sup_fin++;

      if (signe_elt != NULL)
        signe_elt->val[ind_sup] = 1;

    }


  } /* Fin boucle sur les éléments de l'entité supérieure (que l'on fusionne) */

  /* Libération du tableau de recherche */

  ECS_FREE(pos_recherche);
  ECS_FREE(val_recherche);

  /* Compactage de la définition */
  /*-----------------------------*/

  /*
    Le premier élément ne peut pas être fusionné avec un élément précédent,
    la boucle commence donc au deuxième
  */

  cpt_sup_fin       = 1;

  for (ind_sup = 1; ind_sup < table_def->nbr; ind_sup++) {

    if (tab_transf.val[ind_sup] > (ecs_int_t)(cpt_sup_fin - 1)) {

      /* Recopie définition sur partie compactée du tableau */

      for (pos_cpt = table_def->pos[cpt_sup_fin],
             pos_sup = table_def->pos[ind_sup];
           pos_sup < table_def->pos[ind_sup + 1];
           pos_cpt++, pos_sup++)
        table_def->val[pos_cpt - 1] = table_def->val[pos_sup - 1];

      cpt_sup_fin += 1;

      table_def->pos[cpt_sup_fin] = pos_cpt;

    }
  }

  /* Redimensionnement de l'entité compactée */

  table_def->nbr = cpt_sup_fin;
  ECS_REALLOC(table_def->pos, cpt_sup_fin + 1, ecs_size_t);
  ECS_REALLOC(table_def->val,
              table_def->pos[cpt_sup_fin] - 1,
              ecs_int_t);

  ecs_table__pos_en_regle(table_def);

  /* Affectation de la nouvelle liste compactee */
  /*--------------------------------------------*/

  *nbr_elt_cpct = table_def->nbr;

  /* Renvoi du tableau de transformation */
  /*-------------------------------------*/

  return tab_transf;
}

/*----------------------------------------------------------------------------
 *  Fonction qui construit la liste des cellules attachées à une liste
 *  de faces fournie en argument.
 *----------------------------------------------------------------------------*/

ecs_tab_int_t
ecs_table_def__liste_cel_fac(const size_t          nbr_fac,
                             ecs_table_t          *table_def_cel,
                             const ecs_tab_int_t   liste_fac)
{
  size_t      cpt_cel;
  size_t      icel;
  size_t      ifac;
  size_t      iloc;
  size_t      ipos;
  size_t      nbr_cel;
  size_t      nbr_fac_cel;
  ecs_int_t  *indic_cel;
  ecs_int_t  *indic_fac;

  ecs_tab_int_t   liste_cel;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(table_def_cel != NULL);

  ecs_table__regle_en_pos(table_def_cel);

  /* Initialisations */

  cpt_cel = 0;

  nbr_cel = table_def_cel->nbr;

  /* Indicateurs */

  ECS_MALLOC(indic_cel, nbr_cel, ecs_int_t);
  ECS_MALLOC(indic_fac, nbr_fac, ecs_int_t);

  for (icel = 0; icel < nbr_cel; icel++)
    indic_cel[icel] = 0;

  for (ifac = 0; ifac < nbr_fac; ifac++)
    indic_fac[ifac] = 0;

  for (ifac = 0; ifac < liste_fac.nbr; ifac++)
    indic_fac[liste_fac.val[ifac]] = 1;

  /* Première boucle sur les cellules : marquage */

  for (icel = 0; icel < nbr_cel; icel++) {

    nbr_fac_cel
      = table_def_cel->pos[icel + 1]
      - table_def_cel->pos[icel];

    ipos = table_def_cel->pos[icel] - 1;

    for (iloc = 0; iloc < nbr_fac_cel; iloc++) {

      ifac = ECS_ABS(table_def_cel->val[ipos]) - 1;

      ipos++;

      if (indic_fac[ifac] == 1 && indic_cel[icel] == 0) {
        indic_cel[icel] = 1;
        cpt_cel++;
      }

    }
  }

  ECS_FREE(indic_fac);

  /* Seconde boucle sur les cellules : remplissage */

  liste_cel.nbr = cpt_cel;
  ECS_MALLOC(liste_cel.val, liste_cel.nbr, ecs_int_t);

  cpt_cel = 0;

  for (icel = 0; icel < nbr_cel; icel++) {

    if (indic_cel[icel] == 1)
      liste_cel.val[cpt_cel++] = icel;

  }

  ECS_FREE(indic_cel);

  ecs_table__libere_pos(table_def_cel);

  return liste_cel;
}

/*----------------------------------------------------------------------------
 *  Fonction qui remplace les références à des éléments
 *  en des références à d'autres éléments liés aux premiers
 *  par un tableau de renumérotation qui peut être signé.
 *----------------------------------------------------------------------------*/

void
ecs_table_def__remplace_ref(ecs_table_t    *table_def,
                            ecs_tab_int_t  *tab_old_new)
{
  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (table_def != NULL && tab_old_new != NULL) {

    ecs_int_t  num_old, num_new;
    ecs_size_t ielt;
    ecs_size_t ival, ival_deb, ival_fin;
    ecs_size_t cpt_val = 0;
    size_t nbr_elt = table_def->nbr;

    ecs_table__regle_en_pos(table_def);

    ival_deb = 0;

    for (ielt = 0; ielt < nbr_elt; ielt++) {

      ival_fin = table_def->pos[ielt+1] - 1;

      for (ival = ival_deb; ival < ival_fin; ival++) {
        num_new = 0;
        num_old = table_def->val[ival];
        if (num_old > 0)
          num_new = tab_old_new->val[num_old -1];
        else
          num_new = -tab_old_new->val[-num_old -1];
        if (num_new != 0)
          table_def->val[cpt_val++] = num_new;
      }

      ival_deb = table_def->pos[ielt+1] - 1;
      table_def->pos[ielt+1] = cpt_val + 1;
    }

    ECS_REALLOC(table_def->val, cpt_val, ecs_int_t);

    ecs_table__pos_en_regle(table_def);
  }
}

/*----------------------------------------------------------------------------
 *  Fonction qui construit un tableau de booleens conforme a une liste
 *   de sous-elements
 *  Un sous-element est a `true'
 *   s'il intervient dans la definition des elements
 *----------------------------------------------------------------------------*/

void
ecs_table_def__cree_masque(bool          sselt_select[],
                           ecs_table_t  *table_def_elt)
{
  size_t  nbr_elt;
  size_t  num_sselt;
  size_t  pos_inf;
  size_t  pos_sup;

  size_t  ielt;
  size_t  ipos;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  ecs_table__regle_en_pos(table_def_elt);

  nbr_elt = table_def_elt->nbr;

  for (ielt = 0; ielt < nbr_elt; ielt++) {

    pos_inf = table_def_elt->pos[ielt    ] - 1;
    pos_sup = table_def_elt->pos[ielt + 1] - 1;

    for (ipos = pos_inf; ipos < pos_sup; ipos++) {

      num_sselt = ECS_ABS(table_def_elt->val[ipos]) - 1;

      /* If the sub-element has not yet been accounted for, mark it */

      if (sselt_select[num_sselt] == false)
        sselt_select[num_sselt] = true;

    }
  }

  ecs_table__libere_pos(table_def_elt);
}

/*----------------------------------------------------------------------------
 * Suppression des sommets ne participant pas à la connectivité
 *  et mise à jour de la connectivité.
 *----------------------------------------------------------------------------*/

void
ecs_table_def__nettoie_nodal(size_t        *n_vertices,
                             ecs_coord_t  **vtx_coords,
                             ecs_table_t   *table_def_fac,
                             ecs_table_t   *table_def_cel)
{
  size_t  cpt_som;
  size_t  iloc;
  size_t  ipos_new;
  size_t  ipos_old;
  size_t  isom;
  size_t  ival;
  size_t  nbr_som;
  size_t  nbr_val;

  ecs_tab_int_t   tab_som_old_new;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  tab_som_old_new.nbr = 0;
  tab_som_old_new.val = NULL;

  if (vtx_coords == NULL)
    return;

  /* Initialisation du vecteur de renumérotation comme marqueur */

  nbr_som = *n_vertices;

  tab_som_old_new = ecs_tab_int__cree_init(nbr_som, 0);

  /* faces */

  if (table_def_fac != NULL) {

    nbr_val = ecs_table__ret_val_nbr(table_def_fac);

    for (ival = 0; ival < nbr_val; ival++)
      tab_som_old_new.val[table_def_fac->val[ival] - 1] = 1;
  }

  /* cellules */

  if (table_def_cel != NULL) {

    nbr_val = ecs_table__ret_val_nbr(table_def_cel);

    for (ival = 0; ival < nbr_val; ival++)
      tab_som_old_new.val[table_def_cel->val[ival] - 1] = 1;
  }

  /* Préparation de la renumérotation */

  cpt_som = 0;

  for (isom = 0; isom < nbr_som; isom++) {
    if (tab_som_old_new.val[isom] != 0)
      cpt_som++;
  }

  /* Si tous les sommets sont référencés, rien à faire */

  if (cpt_som == nbr_som) {
    tab_som_old_new.nbr = 0;
    ECS_FREE(tab_som_old_new.val);
    return;
  }

  /* Compactage du tableau des sommets
     et transformation du marquer en renumérotation */

  cpt_som = 0;

  for (isom = 0; isom < nbr_som; isom++) {

    if (tab_som_old_new.val[isom] != 0) {

      ipos_new = 3 * (cpt_som);
      ipos_old = 3 * isom;

      for (iloc = 0; iloc < 3; iloc++)
        (*vtx_coords)[ipos_new++] = (*vtx_coords)[ipos_old++];

      tab_som_old_new.val[isom] = cpt_som + 1;

      cpt_som++;

    }
  }

  *n_vertices = cpt_som;
  ECS_REALLOC(*vtx_coords, cpt_som*3, ecs_coord_t);

  /* Mise à jour de la connectivité nodale */

  if (table_def_fac != NULL)
    ecs_table_def__remplace_ref(table_def_fac,
                                &tab_som_old_new);

  if (table_def_cel != NULL)
    ecs_table_def__remplace_ref(table_def_cel,
                                &tab_som_old_new);

  tab_som_old_new.nbr = 0;
  ECS_FREE(tab_som_old_new.val);
}

/*----------------------------------------------------------------------------
 *  Correction si nécessaire de l'orientation des éléments en connectivité
 *   nodale. L'argument liste_cel_err est optionnel.
 *----------------------------------------------------------------------------*/

void
ecs_table_def__orient_nodal(ecs_coord_t     *vtx_coords,
                            ecs_table_t     *table_def_fac,
                            ecs_table_t     *table_def_cel,
                            ecs_tab_int_t   *liste_cel_err,
                            bool             correc_orient)
{
  size_t      ipos_cel;
  size_t      ipos_fac;
  size_t      icel;
  size_t      ifac;
  size_t      nbr_som_loc;
  size_t      nbr_fac;
  size_t      nbr_cel;
  ecs_int_t   typ_elt;
  ecs_int_t   ret_orient;

  ecs_int_t   cpt_orient_correc[ECS_ELT_TYP_FIN];
  ecs_int_t   cpt_orient_erreur[ECS_ELT_TYP_FIN];

  ecs_int_t   typ_cel_base[9] = {ECS_ELT_TYP_NUL,
                                 ECS_ELT_TYP_NUL,
                                 ECS_ELT_TYP_NUL,
                                 ECS_ELT_TYP_NUL,
                                 ECS_ELT_TYP_CEL_TETRA,
                                 ECS_ELT_TYP_CEL_PYRAM,
                                 ECS_ELT_TYP_CEL_PRISM,
                                 ECS_ELT_TYP_NUL,
                                 ECS_ELT_TYP_CEL_HEXA};

  ecs_int_t   cpt_cel_erreur = 0;
  ecs_int_t   cpt_cel_correc = 0;

  int  foam2vtk = 0;  /* tentative adaptation to meshes transformed with
                         foam2vtk; solves most issues but a few remain */
  if (getenv("CS_PREPROCESS_FOAM_2_VTK_SOURCE") != NULL)
    foam2vtk = 1;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (vtx_coords == NULL)
    return;

  if (table_def_fac != NULL)
    ecs_table__regle_en_pos(table_def_fac);

  if (table_def_cel != NULL)
    ecs_table__regle_en_pos(table_def_cel);

  assert(vtx_coords != NULL);

  for (typ_elt = 0; typ_elt < ECS_ELT_TYP_FIN; typ_elt++) {
    cpt_orient_correc[typ_elt] = 0;
    cpt_orient_erreur[typ_elt] = 0;
  }

  printf(_("\n  Element orientation check.\n\n"));

  /* faces */

  if (table_def_fac != NULL) {

    nbr_fac = table_def_fac->nbr;

    for (ifac = 0; ifac < nbr_fac; ifac++) {

      ipos_fac    = table_def_fac->pos[ifac    ] - 1;
      nbr_som_loc = table_def_fac->pos[ifac + 1] - 1 - ipos_fac;

      if (nbr_som_loc == 4) {

        ret_orient
          = _orient_quad(vtx_coords,
                         &(table_def_fac->val[ipos_fac]),
                         correc_orient);

        if (ret_orient < 0)
          cpt_orient_erreur[ECS_ELT_TYP_FAC_QUAD] += 1;
        else if (ret_orient > 0)
          cpt_orient_correc[ECS_ELT_TYP_FAC_QUAD] += 1;

      }
    }
  }

  /* cellules */

  if (table_def_cel != NULL) {

    ecs_tab_int_t face_index, face_marker, edges;

    face_index.nbr = 0;
    face_index.val = NULL;

    face_marker.nbr = 0;
    face_marker.val = NULL;

    edges.nbr = 0;
    edges.val = NULL;

    nbr_cel = table_def_cel->nbr;

    for (icel = 0; icel < nbr_cel; icel++) {

      ipos_cel = table_def_cel->pos[icel] - 1;

      nbr_som_loc = table_def_cel->pos[icel + 1] - 1 - ipos_cel;

      if (nbr_som_loc < 9)
        typ_elt = typ_cel_base[nbr_som_loc];
      else
        typ_elt = ECS_ELT_TYP_CEL_POLY;

      switch (typ_elt) {

      case ECS_ELT_TYP_CEL_TETRA:

        if (foam2vtk)
          ret_orient = _orient_tetra_of(&(table_def_cel->val[ipos_cel]));
        else
          ret_orient = _orient_tetra(vtx_coords,
                                     &(table_def_cel->val[ipos_cel]));
        break;

      case ECS_ELT_TYP_CEL_PYRAM:

        if (foam2vtk)
          ret_orient
            = _orient_pyram_of(&(table_def_cel->val[ipos_cel]));
        else
          ret_orient
            = _orient_pyram(vtx_coords,
                            &(table_def_cel->val[ipos_cel]),
                            correc_orient);
        break;

      case ECS_ELT_TYP_CEL_PRISM:

        if (foam2vtk)
          ret_orient = 0;
        else
          ret_orient
            = _orient_prism(vtx_coords,
                            &(table_def_cel->val[ipos_cel]),
                            correc_orient);
        break;

      case ECS_ELT_TYP_CEL_HEXA:

        if (foam2vtk)
          ret_orient
            = _orient_hexa_of(&(table_def_cel->val[ipos_cel]));
        else
          ret_orient
            = _orient_hexa(vtx_coords,
                           &(table_def_cel->val[ipos_cel]),
                           correc_orient);
        break;

      default: /* ECS_ELT_TYP_CEL_POLY */

        ret_orient = _orient_polyhedron(vtx_coords,
                                        &(table_def_cel->val[ipos_cel]),
                                        (  table_def_cel->pos[icel + 1]
                                         - table_def_cel->pos[icel]),
                                        correc_orient,
                                        &face_index,
                                        &face_marker,
                                        &edges);

      };

      if (ret_orient < 0) {
        cpt_orient_erreur[typ_elt] += 1;
        if (liste_cel_err != NULL) {
          if (liste_cel_err->nbr == 0) {
            liste_cel_err->nbr = nbr_cel;
            ECS_MALLOC(liste_cel_err->val, liste_cel_err->nbr, ecs_int_t);
          }
          liste_cel_err->val[cpt_cel_erreur] = icel;
        }
        cpt_cel_erreur += 1;
      }

      else if (ret_orient > 0) {
        cpt_orient_correc[typ_elt] += 1;
        cpt_cel_correc += 1;
      }
    }

    if (face_index.nbr > 0)
      ECS_FREE(face_index.val);

    if (face_marker.nbr > 0)
      ECS_FREE(face_marker.val);

    if (edges.nbr > 0)
      ECS_FREE(edges.val);
  }

  /* Impression de messages d'avertissement */

  for (typ_elt = 0; typ_elt < ECS_ELT_TYP_FIN; typ_elt++) {
    if (cpt_orient_correc[typ_elt] > 0) {
      ecs_warn();
      printf(_("%d elements of type %s had to be re-oriented\n"),
             (int)(cpt_orient_correc[typ_elt]),
             _(ecs_fic_elt_typ_liste_c[typ_elt].nom));
    }
    if (cpt_orient_erreur[typ_elt] > 0) {
      if (correc_orient == true) {
        ecs_warn();
        printf(_("%d elements of type %s were impossible to re-orient\n"),
               (int)(cpt_orient_erreur[typ_elt]),
               _(ecs_fic_elt_typ_liste_c[typ_elt].nom));
      }
      else {
        ecs_warn();
        printf(_("%d elements of type %s are mis-oriented or highly warped\n"),
               (int)(cpt_orient_erreur[typ_elt]),
               _(ecs_fic_elt_typ_liste_c[typ_elt].nom));
      }
    }
  }

  /* Redimensionnement de la liste des cellules toujours mal orientées */

  if (liste_cel_err != NULL) {
    if (liste_cel_err->nbr > 0) {
      liste_cel_err->nbr = cpt_cel_erreur;
      ECS_REALLOC(liste_cel_err->val, liste_cel_err->nbr, ecs_int_t);
    }
  }

  if (table_def_fac != NULL)
    ecs_table__libere_pos(table_def_fac);

  if (table_def_cel != NULL)
    ecs_table__libere_pos(table_def_cel);
}

/*----------------------------------------------------------------------------
 *  Fusion des sommets confondus d'après la longueur des arêtes des faces.
 * La connectivité des faces est mise à jour.
 *----------------------------------------------------------------------------*/

void
ecs_table_def__nettoie_som_fac(size_t        *n_vertices,
                               ecs_coord_t  **vtx_coords,
                               ecs_table_t   *table_def_fac)
{
  size_t     cpt_fusion;

  size_t     icoo;
  size_t     ifac;

  size_t     ipos_deb;
  size_t     ipos_0;
  size_t     ipos_1;
  size_t     isom;
  size_t     isom_0;
  size_t     isom_1;

  size_t     nbr_fac;
  size_t     nbr_som;
  size_t     nbr_som_fac;

  double      delta_coo;
  double      lng_are_min;
  double      lng_are_2;
  double      lng_min_2;

  ecs_tab_int_t  tab_equiv_som;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (vtx_coords == NULL || table_def_fac == NULL)
    return;

  ecs_table__regle_en_pos(table_def_fac);

  nbr_fac = table_def_fac->nbr;
  nbr_som = *n_vertices;

  lng_are_min = 1.0e-15;
  if (getenv("CS_PREPROCESS_MIN_EDGE_LEN") != NULL) {
    lng_are_min = atof(getenv("CS_PREPROCESS_MIN_EDGE_LEN"));
  }
  lng_min_2 = lng_are_min * ECS_ABS(lng_are_min);

  /* Marquage des sommets */
  /*----------------------*/

  tab_equiv_som = ecs_tab_int__cree_init(nbr_som, -1);

  cpt_fusion = 0;

  /* Boucle sur les faces pour détecter les sommets confondus */
  /*----------------------------------------------------------*/

  for (ifac = 0; ifac < nbr_fac; ifac++) {

    ipos_deb = table_def_fac->pos[ifac] - 1;
    nbr_som_fac = (table_def_fac->pos[ifac + 1] - 1) - ipos_deb;

    for (isom = 0; isom < nbr_som_fac; isom++) {

      isom_0 = table_def_fac->val[ipos_deb + isom] - 1;
      isom_1 = table_def_fac->val[  ipos_deb + ((isom + 1)%nbr_som_fac)] -1;
      ipos_0 = isom_0 * 3;
      ipos_1 = isom_1 * 3;

      lng_are_2 = 0.0;

      for (icoo = 0; icoo < 3; icoo++) {
        delta_coo =   (*vtx_coords)[ipos_1 + icoo]
                    - (*vtx_coords)[ipos_0 + icoo];
        lng_are_2 += delta_coo * delta_coo;
      }

      /* Sommets confondus détectés */

      if (lng_are_2 < lng_min_2) {
        _table_def__maj_equiv_som(isom_0, isom_1, &tab_equiv_som);
        cpt_fusion += 1;
      }

    }
  }

  /* Fusion si nécessaire */

  if (cpt_fusion > 0) {

    ecs_tab_int_t  tab_som_old_new;

    printf(_("\nMesh verification:\n\n"
             "  %d vertices belong to edges of length less than %g\n"
             "  and will be merged; (this tolerance may be modified\n"
             "  using the CS_PREPROCESS_MIN_EDGE_LEN environment variable).\n"),
           (int)cpt_fusion, lng_are_min);

    tab_som_old_new = _table_def__fusion_som(n_vertices,
                                             vtx_coords,
                                             &tab_equiv_som);

    _table_def__maj_fac_som(table_def_fac, &tab_som_old_new);

    ECS_FREE(tab_som_old_new.val);

    printf("\n");
  }

  /* Libération mémoire et retour */

  ECS_FREE(tab_equiv_som.val);

  ecs_table__pos_en_regle(table_def_fac);
}

/*----------------------------------------------------------------------------
 *  Fonction qui supprime les éventuelles faces dégénérées
 *----------------------------------------------------------------------------*/

ecs_tab_int_t
ecs_table_def__nettoie_fac(ecs_table_t  *table_def_fac)
{
  ecs_int_t  nbr_fac_old;
  ecs_int_t  nbr_fac_new;

  ecs_int_t  ind_fac;

  ecs_int_t  cpt_fac;
  ecs_int_t  cpt_val;
  ecs_int_t  ind_val;
  ecs_int_t  ind_val_deb;
  ecs_int_t  ind_val_fin;

  ecs_tab_int_t  tab_fac_old_new;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  tab_fac_old_new.nbr = 0;
  tab_fac_old_new.val = NULL;

  if (table_def_fac == NULL)
    return tab_fac_old_new;

  ecs_table__regle_en_pos(table_def_fac);

  /* Initialisations */

  nbr_fac_old = table_def_fac->nbr;

  /* Boucle sur les faces (renumérotation) */
  /* ------------------------------------- */

  tab_fac_old_new.nbr = nbr_fac_old;
  ECS_MALLOC(tab_fac_old_new.val, tab_fac_old_new.nbr, ecs_int_t);

  cpt_fac = 0;
  cpt_val = 0;

  for (ind_fac = 0; ind_fac < nbr_fac_old; ind_fac++) {

    ind_val_deb = table_def_fac->pos[ind_fac    ] - 1;
    ind_val_fin = table_def_fac->pos[ind_fac + 1] - 1;

    if (ind_val_fin - ind_val_deb > 2) {

      for (ind_val = ind_val_deb; ind_val < ind_val_fin; ind_val++)
        table_def_fac->val[cpt_val++] = table_def_fac->val[ind_val];

      tab_fac_old_new.val[ind_fac] = 1 + cpt_fac++;

      table_def_fac->pos[cpt_fac] = cpt_val + 1;

    }
    else
      tab_fac_old_new.val[ind_fac] = 0;

  }

  nbr_fac_new = cpt_fac;

  /* Si on n'a pas de faces dégénérées, on peut sortir */

  if (nbr_fac_new == nbr_fac_old) {

    tab_fac_old_new.nbr = 0;
    ECS_FREE(tab_fac_old_new.val);

    ecs_table__pos_en_regle(table_def_fac);

    return tab_fac_old_new;
  }

  printf(_("\nMesh verification:\n\n"
           "  Removal of %d degenerate faces:\n"
           "    Initial number of faces           : %10d\n"
           "    Number of faces after processing  : %10d\n"),
         (int)(nbr_fac_old - nbr_fac_new), (int)nbr_fac_old, (int)nbr_fac_new);

  /* On redimensionne le tableau des faces */

  table_def_fac->nbr = nbr_fac_new;
  ECS_REALLOC(table_def_fac->pos, nbr_fac_new + 1, ecs_size_t);
  ECS_REALLOC(table_def_fac->val, cpt_val, ecs_int_t);

  ecs_table__pos_en_regle(table_def_fac);

  return tab_fac_old_new;
}

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie un tableau associant un type à chaque face, sous
 * forme de masque : 0 pour face isolée, 1 ou 2 pour face de bord (1 si
 * cellule avec cette face normale sortante, 2 si cellule avec cette face
 * normale entrante), 1+2 = 3 pour face interne, et 4 ou plus pour tous
 * les autres cas, correspondant à une erreur de connectivité (+4 pour faces
 * voyant au moins deux cellules avec face normale sortante, +8 pour faces
 * voyant au moins deux cellules avec face normale entrante).
 *
 *  Le type de chaque face pourra être modifié ultérieurement en fonction
 * des informations de périodicité.
 *----------------------------------------------------------------------------*/

ecs_tab_int_t
ecs_table_def__typ_fac_cel(ecs_table_t  *table_def_cel,
                           ecs_table_t  *table_def_fac)
{
  size_t      nbr_cel;
  size_t      nbr_fac;
  ecs_int_t   num_fac;
  size_t      icel;
  size_t      ifac;
  size_t      ipos;

  ecs_tab_int_t   typ_fac;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  typ_fac.nbr = 0;
  typ_fac.val = NULL;

  if (table_def_cel == NULL)
    return typ_fac;

  /* Initialisations */

  ecs_table__regle_en_pos(table_def_cel);

  nbr_cel = table_def_cel->nbr;
  nbr_fac = table_def_fac->nbr;

  /* Type de face selon le nombre de cellules voisines : 0 pour face isolée,
     1 ou 2 pour face de bord, 3 pour faces internes, et 4 pour tous les
     autres cas (faces voyant au moins deux cellules sur un même côté) */

  typ_fac.nbr = nbr_fac;

  ECS_MALLOC(typ_fac.val, nbr_fac, ecs_int_t);

  for (ifac = 0; ifac < nbr_fac; ifac++)
    typ_fac.val[ifac] = 0;

  /* Boucle sur les cellules : marquage */
  /*------------------------------------*/

  for (icel = 0; icel < nbr_cel; icel++ ) {

    for (ipos = table_def_cel->pos[icel] - 1;
         ipos < table_def_cel->pos[icel+1] - 1;
         ipos++) {

      num_fac = table_def_cel->val[ipos];

      ifac = ECS_ABS(num_fac) - 1;

      if (num_fac > 0 && (typ_fac.val[ifac] & 1) == 0)
        typ_fac.val[ifac] += 1;

      else if (num_fac < 0 && (typ_fac.val[ifac] & 2) == 0)
        typ_fac.val[ifac] += 2;

      else {
        if (num_fac > 0 && (typ_fac.val[ifac] & 1) == 1)
          typ_fac.val[ifac] = typ_fac.val[ifac] | 4;
        else if (num_fac < 0 && (typ_fac.val[ifac] & 2) == 2)
          typ_fac.val[ifac] = typ_fac.val[ifac] | 8;
      }

    }

  }

  ecs_table__libere_pos(table_def_cel);

  return typ_fac;
}

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie un tableau associant un type à chaque face les
 * numéros des cellules définies par cette face (normale sortante,
 * puis normale entrante). On affecte une valeur 0 lorsqu'il n'y a pas de
 * cellule correspondante directe (la périodicité n'est donc pas prise en
 * compte à ce niveau).
 *
 * On suppose que la cohérence du maillage a déjà été véridifiée et
 * qu'aucune face n'appartient à plus d'une cellule par côté.
 *----------------------------------------------------------------------------*/

ecs_tab_int_t
ecs_table_def__fac_cel(ecs_table_t  *table_def_cel,
                       ecs_table_t  *table_def_fac)
{
  size_t      nbr_cel;
  size_t      nbr_fac;
  ecs_int_t   num_fac;
  size_t      icel;
  size_t      ifac;
  size_t      ipos;

  ecs_tab_int_t   fac_cel;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  fac_cel.nbr = 0;
  fac_cel.val = NULL;

  if (table_def_cel == NULL)
    return fac_cel;

  /* Initialisations */

  ecs_table__regle_en_pos(table_def_cel);

  nbr_cel = table_def_cel->nbr;
  nbr_fac = table_def_fac->nbr;

  /* Allocation et mise à zéro des connectivités */

  fac_cel.nbr = nbr_fac * 2;

  ECS_MALLOC(fac_cel.val, fac_cel.nbr, ecs_int_t);

  for (ipos = 0; ipos < fac_cel.nbr; ipos++)
    fac_cel.val[ipos] = 0;

  /* Boucle sur les cellules : marquage */
  /*------------------------------------*/

  for (icel = 0; icel < nbr_cel; icel++ ) {

    for (ipos = table_def_cel->pos[icel] - 1;
         ipos < table_def_cel->pos[icel+1] - 1;
         ipos++) {

      num_fac = table_def_cel->val[ipos];

      ifac = ECS_ABS(num_fac) - 1;

      if (num_fac > 0) {
        assert(fac_cel.val[ifac*2] == 0);
        fac_cel.val[ifac*2] = icel + 1;
      }
      else {
        assert(fac_cel.val[ifac*2 + 1] == 0);
        fac_cel.val[ifac*2 + 1] = icel + 1;
      }

    }

  }

  ecs_table__libere_pos(table_def_cel);

  return fac_cel;
}

/*----------------------------------------------------------------------------*/

