/*============================================================================
 *  Définition de la fonction de base
 *   de remplissage de la structure de maillage à partir des données lues
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


/*============================================================================
 *                                 Visibilité
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

#include "ecs_def.h"
#include "ecs_mem.h"
#include "ecs_tab.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles des paquetages visibles
 *----------------------------------------------------------------------------*/

#include "ecs_descr.h"
#include "ecs_descr_chaine.h"

#include "ecs_table.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_maillage.h"
#include "ecs_maillage_priv.h"


/*----------------------------------------------------------------------------
 *  Fichier  `include' du  paquetage courant associe au fichier courant
 *----------------------------------------------------------------------------*/

#include "ecs_maillage_pre.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' privés   du  paquetage courant
 *----------------------------------------------------------------------------*/


/*============================================================================
 *                              Fonctions privées
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction qui réalise la transformation
 *   de réferences à un tableau de labels en indices de ce tableau
 *   si le tableau de références a un pas régulier.
 *
 *  Cette fonction retourne :
 *
 *   - `true'  si la table des références a un pas régulier
 *                 (la transformation  a été realisée)
 *   - `false' si la table des références n'a pas un pas régulier
 *                 (aucune transformation n'a été réalisée)
 *----------------------------------------------------------------------------*/

static bool
_ref_reg_en_ind(size_t            nbr_ref,
                size_t            nbr_label,
                ecs_int_t        *val_ref,
                const ecs_int_t  *val_label)
{
  bool        bool_regle;

  size_t      ind;

  ecs_int_t   label_deb;
  ecs_int_t   label_pas;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (nbr_ref == 0)
    return true;

  /* On verifie si les valeurs correspondent a une REGLE */

  bool_regle = true;

  label_deb = val_label[0];

  if (nbr_label > 1) {

    label_pas = val_label[1] - val_label[0];

    for (ind = 2; ind < nbr_label; ind++) {

      if ((val_label[ind] - val_label[ind - 1]) != label_pas) {
        bool_regle = false;
        break;
      }

    }

  }
  else /* if nbr_label == 1 */
    label_pas = 1;


  if (label_pas == 0)
    bool_regle = false;

  if (bool_regle == true) {

    /* La table des labels est régulière */
    /*-----------------------------------*/

    if (label_pas != 1 || label_deb != 1) {
      for (ind = 0; ind < nbr_ref; ind++)
        val_ref[ind] = (val_ref[ind] - label_deb) / label_pas + 1;
    }

    /* else : rien a faire : les références sont déjà égales aux indices */

  }

  return bool_regle;
}

/*----------------------------------------------------------------------------
 *  Fonction qui affiche des informations sur l'entité de maillage
 *   à partir du maillage lu
 *
 *  Les informations affichées sont :
 *  - le nombre d'éléments de l'entité
 *  - les couleurs définies pour cette entité
 *    et pour chaque couleur le nombre d'éléments ayant cette couleur
 *----------------------------------------------------------------------------*/

static void
_ecs_maillage_pre__aff_info_ent(ecs_entmail_t      entite,
                                ecs_int_t          nbr_elt,
                                ecs_int_t          nbr_coul,
                                const ecs_int_t   *coul_val,
                                const ecs_size_t  *nbr_elt_coul)
{
  int    max_coul_val;
  int    max_lng_str_coul;
  int    icoul;
  char   str_num_coul[ECS_STR_SIZE];
  char  *str_coul;

  const char *const ent_msg_info_var_nbr_c[2]
    = {N_("Number of faces"),
       N_("Number of cells")};

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  printf("  ");
  ecs_print_padded_str(_(ent_msg_info_var_nbr_c[entite]), ECS_LNG_AFF_STR);
  printf(" : %*ld\n", ECS_LNG_AFF_ENT, (long)nbr_elt);

  if (nbr_coul != 0) {

    max_coul_val = coul_val[0];

    for (icoul = 1; icoul < nbr_coul; icoul++)
      max_coul_val = ECS_MAX(max_coul_val, coul_val[icoul]);

    /* On affiche au minimum sur 2 tables le numéro de couleur */
    if (max_coul_val < 10)
      max_coul_val = 10;

    sprintf(str_num_coul, "%d", max_coul_val);


    max_lng_str_coul
      = (ecs_int_t)strlen(_("Color"))
      + (ecs_int_t)strlen(str_num_coul)
      + 2; /* Pour ` \0' */

    ECS_MALLOC(str_coul, max_lng_str_coul, char);

    for (icoul = 0; icoul < nbr_coul; icoul++) {

      sprintf(str_coul, "%s %d", _("Color"), (int)(coul_val[icoul]));

      ecs_print_padded_str(NULL, ECS_LNG_AFF_STR - max_lng_str_coul + 2);
      ecs_print_padded_str(str_coul, max_lng_str_coul);
      printf(" : %*ld\n", ECS_LNG_AFF_ENT, (long)(nbr_elt_coul[icoul]));

    }

    ECS_FREE(str_coul);
  }
}

/*----------------------------------------------------------------------------
 *  Fonction qui affiche des informations sur l'entité de maillage
 *   à partir du maillage lu contenant des familles
 *
 *  Les informations affichées sont :
 *  - le nombre d'éléments de l'entité
 *  - pour chaque famille le nombre d'éléments de cette famille
 *----------------------------------------------------------------------------*/

static void
_ecs_maillage_pre__aff_info_fam(ecs_entmail_t  entite,
                                size_t         nbr_elt,
                                const int      fam_val[])
{
  int         max_fam_val;
  int         min_fam_val;
  int         max_lng_str_fam;

  int         ifam;
  size_t      ielt;
  ecs_int_t   nbr_fam_loc;
  ecs_int_t  *nbr_elt_fam;

  char        str_num_fam [ECS_STR_SIZE];
  char       *str_fam;

  const char *const ent_msg_info_var_nbr_c[ECS_N_ENTMAIL]
    = {N_("Number of faces"),
       N_("Number of cells")};

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  printf("  ");
  ecs_print_padded_str(_(ent_msg_info_var_nbr_c[entite]), ECS_LNG_AFF_STR);
  printf(" : %*ld\n", ECS_LNG_AFF_ENT, (long)nbr_elt);

  max_fam_val = fam_val[0];
  min_fam_val = fam_val[0];

  /* Comptage nombre de familles max */

  for (ielt = 1; ielt < nbr_elt; ielt++) {

    if (fam_val[ielt] < min_fam_val)
      min_fam_val = fam_val[ielt];
    else if (fam_val[ielt] > max_fam_val)
      max_fam_val = fam_val[ielt];
  }

  nbr_fam_loc = max_fam_val - min_fam_val + 1;

  ECS_MALLOC(nbr_elt_fam, nbr_fam_loc, ecs_int_t);

  for (ifam = 0; ifam < nbr_fam_loc; ifam++)
    nbr_elt_fam[ifam] = 0;

  /* Comptage nombre de familles par élément */

  for (ielt = 0; ielt < nbr_elt; ielt++)
    nbr_elt_fam[fam_val[ielt] - min_fam_val] += 1;

  /* Affichage */

  if (nbr_fam_loc != 0) {

    /* On affiche au minimum sur 2 tables le numéro de famille */

    sprintf(str_num_fam, "%d", ECS_MAX(10, max_fam_val));

    max_lng_str_fam = (ecs_int_t)(  strlen(_("Family"))
                                  + strlen(str_num_fam)
                                  + 2); /* Pour ` \0' */

    ECS_MALLOC(str_fam, max_lng_str_fam, char);

    for (ifam = 0; ifam < nbr_fam_loc; ifam++) {

      if (nbr_elt_fam[ifam] > 0) {

        sprintf(str_fam, "%s %d", _("Family"), min_fam_val + ifam);

        ecs_print_padded_str(NULL, ECS_LNG_AFF_STR - max_lng_str_fam + 2);
        ecs_print_padded_str(str_fam, max_lng_str_fam);
        printf(" : %*ld\n", ECS_LNG_AFF_ENT, (long)(nbr_elt_fam[ifam]));

      }
    }

    ECS_FREE(str_fam);
  }

  ECS_FREE(nbr_elt_fam);
}

/*============================================================================
 *                             Fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Assign vertex coordinates to a mesh structure.
 *----------------------------------------------------------------------------*/

void
ecs_maillage_pre__cree_som(ecs_maillage_t  *maillage,
                           size_t           nbr_som,
                           ecs_coord_t     *som_val_coord)
{
  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(maillage != NULL);
  assert(maillage->vertex_coords == NULL);
  assert(som_val_coord != NULL);

  /* Set vertex coordinates */

  maillage->n_vertices = nbr_som;
  maillage->vertex_coords = som_val_coord;

  /* Print info */

  printf("  ");
  ecs_print_padded_str(_("Number of vertices"), ECS_LNG_AFF_STR);
  printf(" : %*ld\n", ECS_LNG_AFF_ENT, (long)nbr_som);
}

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie un vecteur de structures "entité de maillage"
 *   pour les éléments (aretes, faces et cellules),
 *   construites à partir de tableaux remplis lors de la lecture des données
 *----------------------------------------------------------------------------*/

void
ecs_maillage_pre__cree_elt(ecs_maillage_t   *maillage,
                           const size_t      nbr_elt_ent[],
                           ecs_size_t       *elt_pos_som_ent[],
                           ecs_int_t        *elt_val_som_ent[],
                           int              *elt_val_famille_ent[],
                           ecs_int_t        *elt_val_couleur_ent[],
                           const ecs_int_t   nbr_coul_ent[],
                           ecs_int_t        *val_coul_ent[],
                           ecs_size_t       *nbr_elt_coul_ent[])
{
  size_t  ielt;
  int     ient;
  int     icoul;

  bool    bool_famille_ent;

  ecs_descr_t  *descr_couleur_tete;
  ecs_descr_t  *descr_couleur;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Boucle de remplissage des entités du maillage */
  /*===============================================*/

  for (ient = 0; ient < ECS_N_ENTMAIL; ient++) {

    /*----------------------------------------------------------*/
    /* S'il y a au moins un élément de cette entité de maillage */
    /*----------------------------------------------------------*/

    if (nbr_elt_ent[ient] != 0) {

      /* Éléments définis en fonction des sommets */
      /*------------------------------------------*/

      assert(elt_val_som_ent[ient] != NULL);

      maillage->table_def[ient] = ecs_table__cree(nbr_elt_ent[ient],
                                                  ECS_PAS_NUL,
                                                  elt_pos_som_ent[ient],
                                                  elt_val_som_ent[ient],
                                                  NULL);

      /* Familles ou couleurs des éléments */
      /*-----------------------------------*/

      bool_famille_ent = false;

      if (elt_val_famille_ent != NULL) {

        if (elt_val_famille_ent[ient] != NULL) {

          bool_famille_ent = true;

          /* Impressions */
          /*-------------*/

          _ecs_maillage_pre__aff_info_fam(ient,
                                          nbr_elt_ent[ient],
                                          elt_val_famille_ent[ient]);

          /* Si les numéros de famille sont tous égaux */
          /*  au numéro de la famille par defaut (`0') */
          /*  on ne crée pas le table famille          */

          ielt = 0;
          while (ielt < nbr_elt_ent[ient] &&
                 elt_val_famille_ent[ient][ielt] == 0)
            ielt++;

          if (ielt != nbr_elt_ent[ient]) {

            /* Tous les numéros de famille ne sont pas égaux   */
            /*  au numéro de la famille par defaut (`0')     : */
            /*  on crée le table famille                       */

            /* Familles des éléments */
            /*-----------------------*/

            maillage->elt_fam[ient]
              = elt_val_famille_ent[ient];

          }
          else
            ECS_FREE(elt_val_famille_ent[ient]);

        }

      }

      if (   bool_famille_ent == false
          && nbr_coul_ent != NULL
          && val_coul_ent != NULL
          && elt_val_couleur_ent != NULL) {

        /* Couleurs des éléments */
        /*-----------------------*/

        /* Création des descripteurs de table "couleur" */

        descr_couleur_tete = NULL;

        for (icoul = 0; icoul < nbr_coul_ent[ient]; icoul++) {

          descr_couleur = ecs_descr__cree(val_coul_ent[ient][icoul],
                                          NULL);

          ecs_descr_chaine__ajoute(&descr_couleur_tete,
                                   descr_couleur);

        }

        if (elt_val_couleur_ent[ient] != NULL) {

          maillage->table_att[ient]
            = ecs_table__cree(nbr_elt_ent[ient],
                              ECS_PAS_UNITE,
                              NULL,
                              elt_val_couleur_ent[ient],
                              descr_couleur_tete);

        }

        /* Impressions */
        /*-------------*/

        _ecs_maillage_pre__aff_info_ent(ient,
                                        nbr_elt_ent[ient],
                                        nbr_coul_ent[ient],
                                        val_coul_ent[ient],
                                        nbr_elt_coul_ent[ient]);

      }

    } /* Fin : s'il y a au moins un élément pour l'entité courante */

    if (val_coul_ent != NULL) {
      if (val_coul_ent[ient] != NULL)
      ECS_FREE(val_coul_ent[ient]);
    }
    if (nbr_elt_coul_ent != NULL) {
      if (nbr_elt_coul_ent[ient] != 0)
        ECS_FREE(nbr_elt_coul_ent[ient]);
    }

  } /* Fin : boucle sur les entités */
}

/*----------------------------------------------------------------------------
 *  Fonction réalisant la transformation de références à des labels en indices
 *----------------------------------------------------------------------------*/

void
ecs_maillage_pre__label_en_indice(size_t      nbr_label,
                                  size_t      nbr_ref,
                                  ecs_int_t  *val_label,
                                  ecs_int_t  *val_ref)
{
  bool       bool_regle;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (nbr_label == 0)
    return;

  assert(val_label != NULL && val_ref != NULL);

  /* Si la table des références est une REGLE */
  /*------------------------------------------*/

  bool_regle = _ref_reg_en_ind(nbr_ref, nbr_label, val_ref, val_label);

  if (bool_regle == false) {

    /* La table des labels n'est pas une REGLE */
    /*-----------------------------------------*/

    /* Il faut rechercher les indices à partir de la table des labels */

    size_t  ind;
    ecs_tab_int_t  tab_label;
    ecs_tab_int_t  tab_ref;

    tab_label.nbr = nbr_label;
    tab_label.val = val_label;

    tab_ref.nbr = nbr_ref;
    tab_ref.val = val_ref;

    tab_ref = ecs_tab_int__ref_en_indice(tab_ref, tab_label, false);

    /* Les indices trouvés par la recherche doivent être incrémentés */
    /*  de `1' car ce sont des positions                             */

    for (ind = 0; ind < nbr_ref; ind++)
      val_ref[ind] += 1;
  }
}

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie l'identificateur de l'entité
 *   auquel appartient le type géométrique donné
 *--------------------------------------------------------------------------- */

ecs_entmail_t
ecs_maillage_pre__ret_typ_geo(int  typ_geo)
{
  ecs_entmail_t entmail = ECS_ENTMAIL_NONE;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Identification de l'entité concernée */
  /*--------------------------------------*/

  switch ((ecs_elt_typ_t)typ_geo) {

  case ECS_ELT_TYP_FAC_TRIA:
  case ECS_ELT_TYP_FAC_QUAD:
  case ECS_ELT_TYP_FAC_POLY:

    entmail = ECS_ENTMAIL_FAC;
    break;

  case ECS_ELT_TYP_CEL_TETRA:
  case ECS_ELT_TYP_CEL_PRISM:
  case ECS_ELT_TYP_CEL_PYRAM:
  case ECS_ELT_TYP_CEL_HEXA:
  case ECS_ELT_TYP_CEL_POLY:

    entmail = ECS_ENTMAIL_CEL;
    break;

  default:

    entmail = ECS_ENTMAIL_NONE;

  }

  return entmail;
}

/*----------------------------------------------------------------------------
 * Fonction qui affiche le nombre d'éléments par entité
 *----------------------------------------------------------------------------*/

void
ecs_maillage_pre__aff_nbr_par_ent(size_t        nbr_som,
                                  const size_t  nbr_elt_ent[],
                                  int           lng_imp)
{

  int max_lng_var_nbr;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  max_lng_var_nbr = ECS_MAX((int)strlen(_("Number of vertices")),
                            (int)strlen(_("Number of faces")) );
  max_lng_var_nbr = ECS_MAX(max_lng_var_nbr,
                            (int)strlen(_("Number of cells")) );

  max_lng_var_nbr = ECS_MAX(max_lng_var_nbr, lng_imp);

  if (nbr_elt_ent[ECS_ENTMAIL_CEL] != 0) {

    printf("    ");
    ecs_print_padded_str(_("Number of cells"), max_lng_var_nbr);
    printf(" : %*lu\n", ECS_LNG_AFF_ENT,
           (unsigned long)(nbr_elt_ent[ECS_ENTMAIL_CEL]));

  }

  if (nbr_elt_ent[ECS_ENTMAIL_FAC] != 0) {

    printf("    ");
    ecs_print_padded_str(_("Number of faces"), max_lng_var_nbr);
    printf(" : %*lu\n", ECS_LNG_AFF_ENT,
           (unsigned long)(nbr_elt_ent[ECS_ENTMAIL_FAC]));

  }

  if (nbr_som != 0) {

    printf("    ");
    ecs_print_padded_str(_("Number of vertices"), max_lng_var_nbr);
    printf(" : %*lu\n", ECS_LNG_AFF_ENT,
           (unsigned long)nbr_som);

  }
}

/*----------------------------------------------------------------------------*/

