/*============================================================================
 *  Définitions des fonctions
 *   associées à la structure `ecs_champ_t' décrivant un champ
 *   et réalisant les sorties pour post-traitement EnSight
 *============================================================================*/

/*
  This file is part of the Code_Saturne Preprocessor, element of the
  Code_Saturne CFD tool.

  Copyright (C) 1999-2010 EDF S.A., France

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
#include <string.h>


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/

#include "ecs_def.h"
#include "ecs_elt_typ_liste.h"
#include "ecs_mem.h"
#include "ecs_file.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage global "Post-Traitement"
 *----------------------------------------------------------------------------*/

#include "ecs_post.h"
#include "ecs_post_ens.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles des paquetages visibles
 *----------------------------------------------------------------------------*/

#include "ecs_famille_chaine.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_champ.h"


/*----------------------------------------------------------------------------
 *  Fichier  `include' du  paquetage courant associé au fichier courant
 *----------------------------------------------------------------------------*/

#include "ecs_post_ens.h"
#include "ecs_champ_post_ens.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' privés   du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_post_ens_priv.h"
#include "ecs_champ_priv.h"


/*============================================================================
 *                              Fonctions privées
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Écriture d'un tableau de coordonnées dans un fichier EnSight Gold.
 *----------------------------------------------------------------------------*/

static void
ecs_loc_post_ens__ecr_coord(const ecs_file_t   *fic,
                            size_t              nbr_som,
                            const ecs_coord_t  *coo_som)
{
  int  icoo;
  size_t  cpt_som, isom, isom_deb;

  float  cnv_val[4096];

  for (icoo = 0; icoo < 3; icoo++) {

    isom_deb = 0;

    while (isom_deb < nbr_som ) {

      for (isom = isom_deb, cpt_som = 0;
           isom < nbr_som && cpt_som < 4096;
           isom++, cpt_som++)
        cnv_val[cpt_som] = coo_som[isom*3 + icoo];

      ecs_file_write(cnv_val, sizeof(float), cpt_som, fic);

      isom_deb = isom;

    }

  }
}

/*----------------------------------------------------------------------------
 *  Écriture d'un tableau d'entiers en tant que réels dans un fichier
 *  EnSight Gold.
 *----------------------------------------------------------------------------*/

static void
ecs_loc_post_ens__ecr_val_champ(const ecs_file_t  *fic,
                                size_t             nbr,
                                const ecs_int_t   *lst,
                                const ecs_int_t   *val)
{
  size_t   ind;
  size_t   ind_deb;
  size_t   nbr_sub;

  const ecs_int_t  *val_ptr;
  float  cnv_val[4096];

  for (ind_deb = 0; ind_deb < nbr; ind_deb += 4096) {

    nbr_sub = 4096;

    if (ind_deb + nbr_sub > nbr)
      nbr_sub = nbr - ind_deb;

    val_ptr = val + ind_deb;

    if (lst != NULL) {

      for (ind = 0; ind < nbr_sub; ind++)
        cnv_val[ind] = val_ptr[lst[ind]];

      ecs_file_write(cnv_val, sizeof(float), nbr_sub, fic);

    }
    else {

      for (ind = 0; ind < nbr_sub; ind++)
        cnv_val[ind] = val_ptr[ind];

      ecs_file_write(cnv_val, sizeof(float), nbr_sub, fic);

    }

  }

}

/*============================================================================
 *                             Fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction qui écrit les connectivités nodales des éléments
 *   selon leur type géometrique
 *
 *  Les élements doivent avoir été triés suivant leur type géometrique
 *----------------------------------------------------------------------------*/

void
ecs_champ_post_ens__ecr_part(const char            *nom_maillage,
                             size_t                 n_vertices,
                             const ecs_coord_t      vertex_coords[],
                             ecs_champ_t           *champ_def,
                             const ecs_tab_int_t   *tab_elt_typ_geo,
                             ecs_post_ens_t        *cas_ens)
{
  size_t     ielt;
  size_t     ival;

  size_t     cpt_elt;
  size_t     cpt_elt_fin;
  size_t     cpt_elt_tot;
  size_t     cpt_buf;
  size_t     nbr_elt;
  size_t     nbr_som;
  int        elt_typ_ref;
  int        marqueur_fin;
  size_t     nbr_elt_typ_geo;
  size_t     nbr_fac_loc;
  size_t     nbr_som_loc;

  ecs_size_t * def_pos_tab;
  ecs_int_t  * def_val_tab;

  ecs_file_t * fic_imp;

  int32_t      connect_buf[4096];

  ecs_post_ens_part_t  * this_part = NULL;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(vertex_coords != NULL);
  assert(champ_def != NULL);

  assert(cas_ens != NULL);

  /* Création de la structure d'information sur la part */
  /*----------------------------------------------------*/

  ECS_REALLOC(cas_ens->tab_part, cas_ens->nbr_part + 1, ecs_post_ens_part_t *);
  cas_ens->nbr_part += 1;

  ECS_MALLOC(this_part, 1, ecs_post_ens_part_t);
  cas_ens->tab_part[cas_ens->nbr_part - 1] = this_part;

  this_part->num_part = cas_ens->nbr_part;
  ECS_MALLOC(this_part->nom_part, strlen(nom_maillage) + 1, char);
  strcpy(this_part->nom_part, nom_maillage);

  this_part->lst_parents = NULL;

  /* Initialisations */
  /*-----------------*/

  fic_imp = ecs_post_ens__ecrit_fic_geo(cas_ens);

  assert(cas_ens->fic_geo != NULL);

  ecs_champ__regle_en_pos(champ_def);

  nbr_som = n_vertices;
  def_pos_tab = champ_def->pos;
  def_val_tab = champ_def->val;

  /* Écriture de l'entête */
  /*----------------------*/

  ecs_post_ens__ecr_chaine(fic_imp, "part");
  ecs_post_ens__ecr_int(fic_imp, this_part->num_part);
  ecs_post_ens__ecr_chaine(fic_imp, this_part->nom_part);

  /* Écriture du nombre de sommets */

  ecs_post_ens__ecr_chaine(fic_imp, "coordinates");
  ecs_post_ens__ecr_int(fic_imp, (int)nbr_som);

  this_part->nbr_som = nbr_som;

  /* Écriture des coordonnées */
  /*--------------------------*/

  /* On est en dimension 3 */

  ecs_loc_post_ens__ecr_coord(fic_imp,
                              nbr_som,
                              vertex_coords);

  /* Écriture des éléments */
  /*------------------------*/

  this_part->nbr_typ_ele = 0;
  this_part->nbr_ele_typ = NULL;
  this_part->nom_typ_ele = NULL;

  cpt_elt = 0;
  nbr_elt = champ_def->nbr;

  cpt_elt_tot = 0;

  elt_typ_ref = -1;

  /* Boucle sur le type d'élément */

  while (cpt_elt < nbr_elt) {

    /* Recherche du prochain type d'élément utilisé */

    elt_typ_ref += 1;

    while (tab_elt_typ_geo->val[elt_typ_ref] == 0)
      elt_typ_ref++;

    cpt_elt_fin = cpt_elt + tab_elt_typ_geo->val[elt_typ_ref];

    nbr_elt_typ_geo = cpt_elt_fin - cpt_elt;

    if (nbr_elt_typ_geo == 0)
      continue;

    /* Écriture du type d'élément */

    switch (elt_typ_ref) {

    case ECS_ELT_TYP_FAC_POLY:
      ecs_post_ens__ecr_chaine(fic_imp, "nsided");
      ecs_post_ens__ecr_int(fic_imp, (int)nbr_elt_typ_geo);
      break;

    case ECS_ELT_TYP_CEL_POLY:
      ecs_post_ens__ecr_chaine(fic_imp, "nfaced");
      ecs_post_ens__ecr_int(fic_imp, (int)nbr_elt_typ_geo);
      break;

    default:
      ecs_post_ens__ecr_chaine(fic_imp,
                               ecs_fic_elt_typ_liste_c[elt_typ_ref].nom);
      ecs_post_ens__ecr_int(fic_imp, (int)nbr_elt_typ_geo);
    }

    this_part->nbr_typ_ele += 1;
    ECS_REALLOC(this_part->nbr_ele_typ, this_part->nbr_typ_ele, ecs_int_t);
    ECS_REALLOC(this_part->nom_typ_ele, this_part->nbr_typ_ele, char *);
    ECS_MALLOC(this_part->nom_typ_ele[this_part->nbr_typ_ele - 1],
               strlen(ecs_fic_elt_typ_liste_c[elt_typ_ref].nom) + 1, char);
    this_part->nbr_ele_typ[this_part->nbr_typ_ele - 1] = nbr_elt_typ_geo;

    switch (elt_typ_ref) {

    case ECS_ELT_TYP_FAC_POLY:
      strcpy(this_part->nom_typ_ele[this_part->nbr_typ_ele - 1],
             "nsided");
      break;

    case ECS_ELT_TYP_CEL_POLY:
      strcpy(this_part->nom_typ_ele[this_part->nbr_typ_ele - 1],
             "nfaced");
      break;

    default:
      strcpy(this_part->nom_typ_ele[this_part->nbr_typ_ele - 1],
             ecs_fic_elt_typ_liste_c[elt_typ_ref].nom);

    }


    /* Écriture des connectivités */
    /*----------------------------*/

    cpt_buf = 0;

    switch(elt_typ_ref) {

    case ECS_ELT_TYP_FAC_POLY: /* Polygones */

      /* Nombre de sommets par élément */

      for (ielt = cpt_elt; ielt < cpt_elt_fin; ielt++) {
        connect_buf[cpt_buf++] = def_pos_tab[ielt + 1] - def_pos_tab[ielt];
        if (cpt_buf == 4096 || ielt + 1 == cpt_elt_fin) {
          ecs_file_write(connect_buf, sizeof(int32_t), cpt_buf, fic_imp);
          cpt_buf = 0;
        }
      }

      /* Connectivité */

      for (ielt = cpt_elt; ielt < cpt_elt_fin; ielt++) {
        cpt_buf = 0;
        for (ival = def_pos_tab[ielt + 1];
             ival > def_pos_tab[ielt    ];
             ival--)
          connect_buf[cpt_buf++] = def_val_tab[ival - 2];
        ecs_file_write(connect_buf, sizeof(int32_t), cpt_buf, fic_imp);
      }

      break;

    case ECS_ELT_TYP_CEL_POLY: /* Polyèdres */

      /* Convention : définition nodale cellule->sommets avec numéros de
         premiers sommets répétés en fin de liste pour marquer la fin
         de chaque face */

      /* Faces par élément */

      for (ielt = cpt_elt; ielt < cpt_elt_fin; ielt++) {
        marqueur_fin = -1;
        nbr_fac_loc = 0;
        for (ival = def_pos_tab[ielt    ] - 1;
             ival < def_pos_tab[ielt + 1] - 1;
             ival++) {
          if (def_val_tab[ival] != marqueur_fin) {
            if (marqueur_fin == -1)
              marqueur_fin = def_val_tab[ival];
          }
          else {
            marqueur_fin = -1;
            nbr_fac_loc += 1;
          }
        }
        connect_buf[cpt_buf++] = nbr_fac_loc;
        if (cpt_buf == 4096 || ielt + 1 == cpt_elt_fin) {
          ecs_file_write(connect_buf, sizeof(int32_t), cpt_buf, fic_imp);
          cpt_buf = 0;
        }
      }

      /* Sommets par face */

      marqueur_fin = -1;
      nbr_som_loc = 0;

      cpt_buf = 0;

      for (ival = def_pos_tab[cpt_elt    ] - 1;
           ival < def_pos_tab[cpt_elt_fin] - 1;
           ival++) {

        if (def_val_tab[ival] != marqueur_fin) {
          nbr_som_loc += 1;
          if (marqueur_fin == -1)
            marqueur_fin = def_val_tab[ival];
        }
        else {
          marqueur_fin = -1;
          connect_buf[cpt_buf++] = nbr_som_loc;
          if (cpt_buf == 4096 || ival + 2 == def_pos_tab[cpt_elt_fin]) {
            ecs_file_write(connect_buf, sizeof(int32_t), cpt_buf, fic_imp);
            cpt_buf = 0;
          }
          nbr_som_loc = 0;
        }

      }

      /* Connectivité */

      marqueur_fin = -1;
      cpt_buf = 0;

      for (ival = def_pos_tab[cpt_elt    ] - 1;
           ival < def_pos_tab[cpt_elt_fin] - 1;
           ival++) {

        if (def_val_tab[ival] != marqueur_fin) {
          connect_buf[cpt_buf++] = def_val_tab[ival];
          if (marqueur_fin == -1)
            marqueur_fin = def_val_tab[ival];
        }
        else {
          marqueur_fin = -1;
          ecs_file_write(connect_buf, sizeof(int32_t), cpt_buf, fic_imp);
          cpt_buf = 0;
        }

      }

      break;

    default: /* Éléments "standard" */

      nbr_som_loc = 0;

      if (cpt_elt < cpt_elt_fin)
        nbr_som_loc = def_pos_tab[cpt_elt + 1] - def_pos_tab[cpt_elt];

      for (ielt = cpt_elt; ielt < cpt_elt_fin; ielt++) {

        for (ival = def_pos_tab[ielt]     - 1;
             ival < def_pos_tab[ielt + 1] - 1;
             ival++)
          connect_buf[cpt_buf++] = def_val_tab[ival];

        if (cpt_buf >= (4096 - nbr_som_loc) || ielt + 1 == cpt_elt_fin) {
          ecs_file_write(connect_buf, sizeof(int32_t), cpt_buf, fic_imp);
          cpt_buf = 0;
        }

      }

      break;

    }


    /* On s'apprête à passer au type d'élément suivant */

    cpt_elt += nbr_elt_typ_geo;

  }

  ecs_file_flush(fic_imp);

  /* Nettoyage avant la sortie */
  /*---------------------------*/

  ecs_champ__libere_pos_tab(champ_def, def_pos_tab);

  /* Fermeture du fichier de géométrie */

  ecs_file_close_stream(cas_ens->fic_geo);
}

/*----------------------------------------------------------------------------
 *  Fonction écrivant le champ à sortir au format Ensight
 *---------------------------------------------------------------------------*/

void
ecs_champ_post_ens__ecr_val(const ecs_tab_int_t  *tab_val,
                            const char           *nom_maillage,
                            const char           *nom_champ,
                            ecs_post_ens_t       *cas_ens)
{
  size_t     ielt;
  ecs_int_t  ind_part;

  ecs_int_t  ielt_max;
  ecs_int_t  ind_typ_ele;

  ecs_post_ens_part_t * part_loc;

  ecs_file_t        * fic_champ;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(tab_val != NULL);

  /* Recherche du `part' correspondant */
  /*-----------------------------------*/

  for (ind_part = 0;
       (   ind_part < cas_ens->nbr_part
        && strcmp((cas_ens->tab_part[ind_part])->nom_part,
                  nom_maillage) != 0);
       ind_part++);

  if (ind_part >= cas_ens->nbr_part)
    return;

  fic_champ = ecs_post_ens__ecrit_fic_var(cas_ens, nom_champ);

  while (   ind_part < cas_ens->nbr_part
         && strcmp((cas_ens->tab_part[ind_part])->nom_part,
                   nom_maillage) == 0) {

    part_loc = cas_ens->tab_part[ind_part];
    ind_part++;

    /* Écriture de l'entête */
    /*----------------------*/

    ecs_post_ens__ecr_chaine(fic_champ, "part");

    ecs_post_ens__ecr_int(fic_champ, part_loc->num_part);

    /* Initialisations */
    /*-----------------*/

    ielt_max = -1;

    /* Boucle sur le type d'élément */
    /*------------------------------*/

    ielt = 0;

    for (ind_typ_ele = 0;
         ind_typ_ele < part_loc->nbr_typ_ele;
         ind_typ_ele++) {

      if (   strcmp(part_loc->nom_typ_ele[ind_typ_ele], "nsided") != 0
          && strcmp(part_loc->nom_typ_ele[ind_typ_ele], "nfaced") != 0) {

        ecs_post_ens__ecr_chaine(fic_champ, part_loc->nom_typ_ele[ind_typ_ele]);

        /* Écriture des valeurs */

        ecs_loc_post_ens__ecr_val_champ(fic_champ,
                                        part_loc->nbr_ele_typ[ind_typ_ele],
                                        part_loc->lst_parents,
                                        tab_val->val + ielt);

      }

      ielt += part_loc->nbr_ele_typ[ind_typ_ele];

    }

  } /* Fin de la boucle sur les parts associés au maillage */

  ecs_file_free(fic_champ);
}

/*----------------------------------------------------------------------------*/

