/*============================================================================
 *  Definitions des fonctions
 *   associees a la structure `ecs_maillage_t' decrivant un maillage
 *   et realisant les sorties pour post-traitement
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

#include "cs_config.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie standard C
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <string.h>


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/

#include "ecs_def.h"
#include "ecs_fic.h"
#include "ecs_mem.h"
#include "ecs_tab.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage global "Post-Traitement"
 *----------------------------------------------------------------------------*/

#include "ecs_post_ens.h"
#include "ecs_post.h"

#if defined(HAVE_CGNS)
#include "ecs_post_cgns.h"
#endif /* HAVE_CGNS */

#if defined(HAVE_MED)
#include "ecs_post_med.h"
#endif /* HAVE_MED */


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles des paquetages visibles
 *----------------------------------------------------------------------------*/

#include "ecs_table.h"
#include "ecs_table_att.h"
#include "ecs_table_def.h"
#include "ecs_table_post.h"
#include "ecs_famille.h"

#include "ecs_maillage.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage courant
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichier  `include' du  paquetage courant associe au fichier courant
 *----------------------------------------------------------------------------*/

#include "ecs_maillage_post.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' prives   du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_maillage_priv.h"


/*============================================================================
 *                              Fonctions privees
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction écrivant les connectivités (nodales) des élements
 *----------------------------------------------------------------------------*/

static void
_maillage_post__ecr_elt(const char            *nom_maillage,
                        const ecs_maillage_t  *maillage,
                        ecs_post_type_t        type_post,
                        ecs_post_t            *cas_post)
{
  int             dim_entite_max = 2;
  ecs_table_t    *table_def_sup = NULL;
  ecs_table_t    *table_def_inf = NULL;
  int            *elt_fam_sup = NULL;
  int            *elt_fam_inf = NULL;
  ecs_famille_t  *famille_sup = NULL;
  ecs_famille_t  *famille_inf = NULL;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (cas_post == NULL)
    return;

  assert(maillage != NULL);

  /* Connectivité et familles inférieures pour les cellules */

  if (maillage->table_def[ECS_ENTMAIL_CEL] != NULL) {

    dim_entite_max = 3;

    table_def_sup = maillage->table_def[ECS_ENTMAIL_CEL];
    elt_fam_sup = maillage->elt_fam[ECS_ENTMAIL_CEL];
    famille_sup = maillage->famille[ECS_ENTMAIL_CEL];

    table_def_inf = maillage->table_def[ECS_ENTMAIL_FAC];

    if (maillage->elt_fam[ECS_ENTMAIL_FAC] != NULL) {
      elt_fam_inf = maillage->elt_fam[ECS_ENTMAIL_FAC];
      famille_inf = maillage->famille[ECS_ENTMAIL_FAC];
    }

  }
  else {

    dim_entite_max = 2;

    table_def_sup = maillage->table_def[ECS_ENTMAIL_FAC];

    if (maillage->elt_fam[ECS_ENTMAIL_FAC] != NULL) {
      elt_fam_sup = maillage->elt_fam[ECS_ENTMAIL_FAC];
      famille_sup = maillage->famille[ECS_ENTMAIL_FAC];
    }

  }

  /* Écriture du maillage */

  ecs_table_post__ecr_elt(nom_maillage,
                          dim_entite_max,
                          maillage->n_vertices,
                          maillage->vertex_coords,
                          table_def_sup,
                          elt_fam_sup,
                          table_def_inf,
                          elt_fam_inf,
                          famille_sup,
                          famille_inf,
                          type_post,
                          cas_post);
}

/*============================================================================
 *                             Fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *    Fonction qui écrit le maillage sur fichier pour post-traitement.
 *
 *    Les éléments de l'entité à écrire doivent être sous forme de
 *  connectivité nodale.
 *----------------------------------------------------------------------------*/

void
ecs_maillage_post__ecr(const char       *nom_maillage,
                       ecs_maillage_t   *maillage,
                       ecs_post_type_t   type_post,
                       ecs_post_t       *cas_post)
{
  bool  bool_ecrit_maillage = false;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(maillage != NULL);

  /* On vérifie si l'on a besoin d'écrire le maillage */
  /*--------------------------------------------------*/

  if (cas_post->opt_ens[type_post] == true)
    bool_ecrit_maillage = true;

#if defined(HAVE_CGNS)
  if (cas_post->opt_cgns[type_post] == true)
    bool_ecrit_maillage = true;
#endif

#if defined(HAVE_MED)
  if (cas_post->opt_med[type_post] == true)
    bool_ecrit_maillage = true;
#endif

  if (bool_ecrit_maillage == false)
    return;

  /* Ajout du cas EnSight (si nécessaire) */
  /*--------------------------------------*/

  if (cas_post->opt_ens[type_post] == true) {

    printf(_("\n\n"
             "EnSight output of mesh: %s\n"
             "-----------------------\n"), nom_maillage);

    if (cas_post->cas_ens == NULL)
      cas_post->cas_ens = ecs_post_ens__cree_cas(cas_post->nom_cas);

  }

#if defined(HAVE_CGNS)

  /* Ajout du cas CGNS (si nécessaire) */
  /*-----------------------------------*/

  if (cas_post->opt_cgns[type_post] == true) {

    printf(_("\n\n"
             "CGNS file output of mesh: %s\n"
             "-------------------------\n"), nom_maillage);


    if (cas_post->cas_cgns == NULL)
      cas_post->cas_cgns = ecs_post_cgns__cree_cas(cas_post->nom_cas);
  }

#endif /* HAVE_CGNS */

#if defined(HAVE_MED)

  /* Ajout du cas MED (si nécessaire) */
  /*----------------------------------*/

  if (cas_post->opt_med[type_post] == true) {

    printf(_("\n\n"
             "MED file output of mesh: %s\n"
             "------------------------\n"), nom_maillage);

    if (cas_post->cas_med == NULL)
      cas_post->cas_med = ecs_post_med__cree_cas(cas_post->nom_cas);

  }

#endif /* HAVE_MED */

  /* Appel a la fonction réalisant l'ecriture des élements */
  /*-------------------------------------------------------*/

  if (   maillage->table_def[ECS_ENTMAIL_FAC] != NULL
      || maillage->table_def[ECS_ENTMAIL_CEL] != NULL)
    _maillage_post__ecr_elt(nom_maillage,
                            maillage,
                            type_post,
                            cas_post);
}

/*----------------------------------------------------------------------------
 *  Écriture du maillage correspondant à une liste de faces sur fichier
 *  pour post-traitement.
 *
 *  Cette fonction crée une coupe correspondant à la liste de faces donnée
 *  (ce qui provoque automatiquement son post-traitement), puis la détruit.
 *  Le nom utilisé pour cette sortie ne sera donc plus disponible pour
 *  d'autres coupes.
 *
 *  Le maillage principal doit être en connectivité descendante.
 *----------------------------------------------------------------------------*/

void
ecs_maillage_post__ecr_fac_liste(const char           *nom_liste,
                                 ecs_maillage_t       *maillage,
                                 const ecs_tab_int_t   liste_fac,
                                 ecs_post_type_t       type_post,
                                 ecs_post_t           *cas_post)
{
  ecs_maillage_t  * maillage_extrait;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(maillage != NULL);

  /* Initialisations */
  /*-----------------*/

  maillage_extrait = NULL;

  /* Si aucun post-traitement, rien à faire */

  if (cas_post == NULL)
    return;

  /* Extraction des faces sélectionnées */
  /*------------------------------------*/

  maillage_extrait = ecs_maillage__extrait(maillage,
                                           ECS_ENTMAIL_FAC,
                                           &liste_fac);

  maillage_extrait->typ_connect = maillage->typ_connect;


  /* Tri selon le type des faces */
  /*-----------------------------*/

  printf(_("\nCreating mesh for output: %s"
           "\n-------------------------\n\n"),
         nom_liste);

  ecs_maillage__trie_typ_geo(maillage_extrait);

  /* Sortie de la géométrie pour post-traitement */

  ecs_maillage_post__ecr(nom_liste,
                         maillage_extrait,
                         type_post,
                         cas_post);

  /* Suppression du maillage extrait */

  ecs_maillage__detruit(&maillage_extrait);
}

/*----------------------------------------------------------------------------*/

