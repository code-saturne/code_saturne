/*============================================================================
 *  Définitions des fonctions
 *   associées à la structure `ecs_table_t' décrivant une table
 *   et réalisant les sorties pour post-traitement
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

#include "cs_config.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie standard C
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <string.h> /* strlen() */


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage global "Utilitaire"
 *---------------------------------------------------------------------------*/

#include "ecs_def.h"
#include "ecs_file.h"
#include "ecs_mem.h"
#include "ecs_tab.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage global "Post-Traitement"
 *---------------------------------------------------------------------------*/

#include "ecs_post_ens.h"
#if defined(HAVE_CGNS)
#include "ecs_post_cgns.h"
#endif /* HAVE_CGNS */
#if defined(HAVE_MED)
#include "ecs_post_med.h"
#endif /* HAVE_MED */


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles des paquetages visibles
 *---------------------------------------------------------------------------*/

#include "ecs_table_post_ens.h"

#if defined(HAVE_CGNS)
#include "ecs_table_post_cgns.h"
#endif /* HAVE_CGNS */
#if defined(HAVE_MED)
#include "ecs_table_post_med.h"
#endif /* HAVE_MED */


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage courant
 *---------------------------------------------------------------------------*/

#include "ecs_table.h"


/*----------------------------------------------------------------------------
 *  Fichier  `include' du  paquetage courant associe au fichier courant
 *---------------------------------------------------------------------------*/

#include "ecs_table_post.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' prives   du  paquetage courant
 *---------------------------------------------------------------------------*/

#include "ecs_table_priv.h"
#include "ecs_post_ens_priv.h"

#if defined(HAVE_CGNS)
#include "ecs_post_cgns_priv.h"
#endif

#if defined(HAVE_MED)
#include "ecs_med_priv.h"
#endif


/*============================================================================
 *                       Prototypes de fonctions privées
 *============================================================================*/

/*============================================================================
 *                             Fonctions privées
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction qui compte le nombre d'éléments de chaque type géométrique.
 *  On renvoie un tableau donnant le nombre d'éléments de chaque type
 *  décrit dans "ecs_elt_typ_glob.h"
 *
 *  Les éléments doivent avoir ete triés suivant leur type géometrique.
 *---------------------------------------------------------------------------*/

static ecs_tab_int_t
ecs_loc_table_post__cpt_elt_typ(const ecs_table_t  *table_def,
                                int                 dim_entite)
{
  size_t       cpt_elt;
  size_t       cpt_elt_typ_geo;
  ecs_int_t    elt_typ_ref;
  ecs_int_t    elt_typ_ref_prec;
  size_t       ielt;
  size_t       nbr_val_ref;

  size_t       nbr_elt;
  ecs_size_t * elt_def_pos;
  ecs_int_t  * elt_typ_geo;

  ecs_tab_int_t  nbr_elt_typ_geo;

  const ecs_elt_typ_t  * typ_geo_base;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert (table_def != NULL);
  assert(dim_entite >=1 && dim_entite <= 3);

  /* Initialisations */

  cpt_elt = 0;
  nbr_elt = 0;

  nbr_elt_typ_geo = ecs_tab_int__cree_init(ECS_ELT_TYP_FIN, 0);

  elt_typ_ref = -1;
  elt_typ_ref_prec = -1;

  elt_def_pos = NULL;
  elt_typ_geo = NULL;

  /* Détermination de base du type d'élément "classique" */

  typ_geo_base = ecs_glob_typ_elt[dim_entite - 2];

  /* On utilise les tailles des éléments */
  /*-------------------------------------*/

  nbr_elt     = table_def->nbr;
  elt_def_pos = table_def->pos;

  if (elt_def_pos == NULL) {

    if (table_def->pas < 9)
      elt_typ_ref = typ_geo_base[table_def->pas];
    else if (dim_entite == 2)
      elt_typ_ref = ECS_ELT_TYP_FAC_POLY;
    else
      elt_typ_ref = ECS_ELT_TYP_CEL_POLY;

    nbr_elt_typ_geo.val[elt_typ_ref] = nbr_elt;

  }
  else {

    /* Comptage du nombre d'éléments de chaque type */

    elt_typ_ref_prec = ECS_ELT_TYP_NUL;

    while (cpt_elt < nbr_elt) {

      cpt_elt_typ_geo = 0;

      for (ielt = cpt_elt; ielt < nbr_elt; ielt++) {

        nbr_val_ref = elt_def_pos[ielt + 1] - elt_def_pos[ielt];

        if (nbr_val_ref < 9)
          elt_typ_ref = typ_geo_base[nbr_val_ref];
        else if (dim_entite == 2)
          elt_typ_ref = ECS_ELT_TYP_FAC_POLY;
        else
          elt_typ_ref = ECS_ELT_TYP_CEL_POLY;

        if (elt_typ_ref == elt_typ_ref_prec)
          cpt_elt_typ_geo++;

        else
          break;

      }

      cpt_elt += cpt_elt_typ_geo;

      /* Si les éléments sont bien triés, on ne doit rencontrer chaque
         type que dans un sous-ensemble contigu. */

      if (elt_typ_ref < elt_typ_ref_prec)
        ecs_error(__FILE__, __LINE__, 0,
                  _("Elements of the mesh to post-process have\n"
                    "not been sorted by type."));

      /* Mise à jour du tableau à renvoyer */

      if (cpt_elt_typ_geo > 0)
        nbr_elt_typ_geo.val[elt_typ_ref_prec] = cpt_elt_typ_geo;

      elt_typ_ref_prec = elt_typ_ref;

    }

  }

  /* Retour */

  return nbr_elt_typ_geo;
}

/*============================================================================
 *                             Fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction ecrivant les elements d'une table donne pour le post traitement
 *
 *  Les elements doivent avoir ete tries suivant leur type geometrique
 *---------------------------------------------------------------------------*/

void
ecs_table_post__ecr_elt(const char            *nom_maillage,
                        int                    dim_entite_max,
                        size_t                 n_vertices,
                        ecs_coord_t            vertex_coords[],
                        ecs_table_t           *table_def,
                        const int              elt_fam[],
                        ecs_table_t           *table_def_inf,
                        const int              elt_fam_inf[],
                        const ecs_famille_t   *famille_elt,
                        const ecs_famille_t   *famille_inf,
                        ecs_post_type_t        type_post,
                        ecs_post_t            *cas_post)
{
  ecs_tab_int_t     tab_elt_typ_geo;
  ecs_tab_int_t     tab_elt_typ_geo_inf;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (cas_post == NULL)
    return;

  assert(vertex_coords != NULL);
  assert(table_def != NULL);

  tab_elt_typ_geo = ecs_loc_table_post__cpt_elt_typ(table_def,
                                                    dim_entite_max);

  /* Les elements ayant deja ete tries suivant leur type geometrique ... */
  /* ------------------------------------------------------------------- */

  if (   cas_post->cas_ens != NULL
      && cas_post->opt_ens[type_post] == true) {

    ecs_table_post_ens__ecr_part(nom_maillage,
                                 n_vertices,
                                 vertex_coords,
                                 table_def,
                                 &tab_elt_typ_geo,
                                 cas_post->cas_ens);

  }

#if defined(HAVE_CGNS)

  if (   cas_post->cas_cgns != NULL
      && cas_post->opt_cgns[type_post] == true) {

    ecs_post_cgns__ajoute_maillage(nom_maillage,
                                   dim_entite_max,
                                   cas_post->cas_cgns);

    ecs_table_post_cgns__ecr_connect(nom_maillage,
                                     n_vertices,
                                     vertex_coords,
                                     table_def,
                                     elt_fam,
                                     &tab_elt_typ_geo,
                                     cas_post->cas_cgns);

    ecs_post_cgns__ferme_cas(cas_post->cas_cgns);

  }

#endif /* HAVE_CGNS */

#if defined(HAVE_MED)

  if (   cas_post->cas_med != NULL
      && cas_post->opt_med[type_post] == true) {

    ecs_post_med__ajoute_maillage(nom_maillage,
                                  dim_entite_max,
                                  cas_post->cas_med);

    /* Écriture des familles */

    ecs_table_post_med__ecr_famille(nom_maillage,
                                    famille_elt,
                                    famille_inf,
                                    cas_post->cas_med);

    /* Écriture du maillage principal */

    ecs_table_post_med__ecr_som(nom_maillage,
                                n_vertices,
                                vertex_coords,
                                cas_post->cas_med);

    ecs_table_post_med__ecr_elt(nom_maillage,
                                table_def,
                                elt_fam,
                                &tab_elt_typ_geo,
                                cas_post->cas_med);


    /* Écriture des faces (et de leurs familles) */

    if (table_def_inf != NULL) {

      tab_elt_typ_geo_inf = ecs_loc_table_post__cpt_elt_typ(table_def_inf,
                                                            dim_entite_max-1);

      ecs_table_post_med__ecr_elt(nom_maillage,
                                  table_def_inf,
                                  elt_fam_inf,
                                  &tab_elt_typ_geo_inf,
                                  cas_post->cas_med);

      tab_elt_typ_geo_inf.nbr = 0;
      ECS_FREE(tab_elt_typ_geo_inf.val);
    }

    /* Mise à jour de la structure maillage_med pour y indiquer le nombre
       d'éléments principaux de chaque type MED */

    ecs_table_post_med__cpt_elt_typ(&tab_elt_typ_geo,
                                    nom_maillage,
                                    cas_post->cas_med);

  }

#endif /* HAVE_MED */

  tab_elt_typ_geo.nbr = 0;
  ECS_FREE(tab_elt_typ_geo.val);
}

/*----------------------------------------------------------------------------
 *  Fonction ecrivant les valeurs d'une table donne pour le post traitement
 *---------------------------------------------------------------------------*/

void
ecs_table_post__ecr_val(const ecs_tab_int_t  *tab_val,
                        const char           *nom_maillage,
                        const char           *nom_table,
                        ecs_post_t           *cas_post)
{
  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(tab_val != NULL);

  if (cas_post->cas_ens != NULL) {

    ecs_table_post_ens__ecr_val(tab_val,
                                nom_maillage,
                                nom_table,
                                cas_post->cas_ens);

  }

#if defined(HAVE_CGNS)

  if (cas_post->cas_cgns != NULL) {

    /* Sortie de la table sur fichier CGNS */

    ecs_table_post_cgns__ecr_val(tab_val,
                                 nom_maillage,
                                 nom_table,
                                 cas_post->cas_cgns);
  }

#endif /* HAVE_CGNS */

#if defined(HAVE_MED)

  if (cas_post->cas_med != NULL) {

    ecs_int_t    ind;

    med_type_champ type_champ_med = MED_INT32;
    med_err        ret_med        = 0;

    /* Création du champ si première instance */

    if (ecs_post_med__test_champ_liste(nom_table,
                                       cas_post->cas_med) == false) {

      char   nom_champ_med[MED_TAILLE_NOM + 1];
      char  *nom_comp_med;
      char  *uni_comp_med;

      /* Nom du champ */

      strncpy(nom_champ_med, nom_table, MED_TAILLE_NOM);
      nom_champ_med[MED_TAILLE_NOM] = '\0';

      if (ecs_post_med__test_champ_liste(nom_table,
                                         cas_post->cas_med) == false) {

        /* Nom et unite */

        ECS_MALLOC(nom_comp_med, MED_TAILLE_PNOM + 1, char);
        ECS_MALLOC(uni_comp_med, MED_TAILLE_PNOM + 1, char);
        for (ind = 0; ind < (ecs_int_t)(MED_TAILLE_PNOM); ind++)
          nom_comp_med[ind] = ' ', uni_comp_med[ind] = ' ';
        nom_comp_med[MED_TAILLE_PNOM] = '\0';
        uni_comp_med[MED_TAILLE_PNOM] = '\0';

        ecs_post_med__ajoute_champ_liste(nom_table, cas_post->cas_med);

        if (sizeof(med_int) == 8)
          type_champ_med = MED_INT64;

        ret_med = MEDchampCr(cas_post->cas_med->fid,
                             nom_champ_med,
                             type_champ_med,
                             nom_comp_med,
                             uni_comp_med,
                             1);

        if (ret_med != 0)
          ecs_error(__FILE__, __LINE__, 0,
                    _("MED: error writing file \"%s\".\n"
                      "Error creating field \"%s\"\n"),
                    cas_post->cas_med->nom_fic, nom_champ_med);

        ECS_FREE(nom_comp_med);
        ECS_FREE(uni_comp_med);

      }
    }

    /* Ecriture des valeurs de la table */

    ecs_table_post_med__ecr_val(tab_val,
                                nom_maillage,
                                nom_table,
                                cas_post->cas_med);
  }

#endif /* HAVE_MED */
}

/*----------------------------------------------------------------------------*/

