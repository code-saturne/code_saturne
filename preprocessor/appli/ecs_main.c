/*============================================================================
 *  Programme principal de l'Enveloppe du Code_Saturne
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
 *                                 Visibilité
 *============================================================================*/

#include "cs_config.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie standard C
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <locale.h>
#include <stdlib.h>
#include <string.h>


/*----------------------------------------------------------------------------
 *  Fichiers `include' système
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/

#include "ecs_comm.h"
#include "ecs_def.h"
#include "ecs_file.h"
#include "ecs_mem.h"
#include "ecs_mem_usage.h"
#include "ecs_tab.h"
#include "ecs_timer.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage global "Pre-Post-Traitement"
 *----------------------------------------------------------------------------*/

#include "ecs_post.h"
#include "ecs_pre.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles des paquetages visibles
 *----------------------------------------------------------------------------*/

#include "ecs_maillage.h"
#include "ecs_maillage_post.h"
#include "ecs_maillage_ncs.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_cmd.h"


/*----------------------------------------------------------------------------
 *  Fichier `include' du paquetage courant associé au fichier courant
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' privés du  paquetage courant
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Variables globales statiques
 *----------------------------------------------------------------------------*/


/*============================================================================
 *                              Fonctions privées
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction affichant la taille théorique d'un maillage
 *----------------------------------------------------------------------------*/

static void
_aff_taille_maillage(const ecs_maillage_t *maillage)
{
  size_t       pr_size;
  double       f_size;
  int          nb_div;
  static char  unit_prefix[] = {' ', 'K', 'M', 'G'};

  f_size = ecs_maillage__ret_taille(maillage);

  for (nb_div = 0; f_size > 1024. && nb_div < 3; nb_div++)
    f_size /= 1024.;

  printf(_("  Theoretical mesh size:       %15.3f %cb\n"),
         f_size, unit_prefix[nb_div]);

  f_size = ecs_mem_size_current();

  for (nb_div = 1; f_size > 1024.0 && nb_div < 3; nb_div++)
    f_size /= 1024.;

  printf(_("  Theoretical current memory:  %15.3f %cb\n"),
         f_size, unit_prefix[nb_div]);

  f_size = ecs_mem_size_max();

  for (nb_div = 1; f_size > 1024.0 && nb_div < 3; nb_div++)
    f_size /= 1024.;

  printf(_("  Theoretical peak memory:     %15.3f %cb\n"),
         f_size, unit_prefix[nb_div]);

  pr_size = ecs_mem_usage_max_pr_size();

  if (pr_size > 0) {

    f_size = pr_size;

    for (nb_div = 1; f_size > 1024.0 && nb_div < 3; nb_div++)
      f_size /= 1024.;

    printf(_("  Total memory used:           %15.3f %cb\n"),
           f_size, unit_prefix[nb_div]);

  }
}

/*----------------------------------------------------------------------------
 *  Fonction qui lit les maillages sur fichiers
 *
 *  La fonction renvoie le maillage concaténé
 *----------------------------------------------------------------------------*/

static ecs_maillage_t *
_lit_maillage(const ecs_cmd_t *cmd)
{
  int              ific;

  ecs_maillage_t  *maillage = NULL;
  ecs_maillage_t  *maillage_lu = NULL;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for (ific = 0; ific < cmd->n_num_maillage; ific++) {

    /* Chargement de la structure */

    maillage_lu = ecs_pre__lit_maillage(cmd->fic_maillage,
                                        cmd->fmt_maillage,
                                        cmd->num_maillage[ific],
                                        cmd->grp_cel_section,
                                        cmd->grp_cel_zone,
                                        cmd->grp_fac_section,
                                        cmd->grp_fac_zone);

    if (ific == 0)
      maillage = maillage_lu;

    else
      ecs_maillage__concatene_nodal(maillage, maillage_lu);

  } /* Fin : boucle sur les fichiers de maillage à lire */

  if (cmd->n_num_maillage == 1)
    printf (_("\nDone reading mesh"
              "\n-----------------\n"));
  else
    printf (_("\nDone reading and concatenating meshes"
              "\n-------------------------------------\n"));

  if (cmd->nbr_dump > 0)
    ecs_maillage__imprime(maillage, cmd->nbr_dump);

  _aff_taille_maillage(maillage);

  /* Retour du maillage                                        */
  /*  construit à partir de la concaténation des maillages lus */

  return maillage;
}

/*----------------------------------------------------------------------------
 *  Fonction qui prépare un cas de post traitement si nécessaire
 *----------------------------------------------------------------------------*/

static ecs_post_t *
_init_post(const ecs_cmd_t *cmd)
{
  ecs_post_t     *cas_post;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  cas_post = ecs_post__cree_cas(cmd->nom_cas);

  if (!strcmp(cmd->post_err, "ens"))
    cas_post->opt_ens[ECS_POST_TYPE_ERREUR] = true;
  if (!strcmp(cmd->post_vol, "ens"))
    cas_post->opt_ens[ECS_POST_TYPE_VOLUME] = true;

#if defined(HAVE_CGNS)
  if (!strcmp(cmd->post_err, "cgns"))
    cas_post->opt_cgns[ECS_POST_TYPE_ERREUR] = true;
  if (!strcmp(cmd->post_vol, "cgns"))
    cas_post->opt_cgns[ECS_POST_TYPE_VOLUME] = true;
#endif

#if defined(HAVE_MED)
  if (!strcmp(cmd->post_err, "med"))
    cas_post->opt_med[ECS_POST_TYPE_ERREUR] = true;
  if (!strcmp(cmd->post_vol, "med"))
    cas_post->opt_med[ECS_POST_TYPE_VOLUME] = true;
#endif

  return cas_post;
}

/*----------------------------------------------------------------------------
 *  Impression de cellules (ou faces pour une maillage de peau uniquement)
 *   correspondant à un critère de sélection donné
 *----------------------------------------------------------------------------*/

static void
_post_ele_liste(ecs_maillage_t       *maillage,
                const ecs_tab_int_t   liste_filtre,
                const char           *nom,
                ecs_post_type_t       type_post,
                ecs_post_t           *cas_post)
{
  ecs_maillage_t  *maillage_post;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  maillage_post
    = ecs_maillage__extrait(maillage,
                            (   ecs_maillage__ret_entmail_max(maillage)
                             == ECS_ENTMAIL_CEL ?
                                ECS_ENTMAIL_CEL : ECS_ENTMAIL_FAC),
                            &liste_filtre);

  ecs_maillage_post__ecr(nom,
                         maillage_post,
                         type_post,
                         cas_post);

  ecs_maillage__detruit(&maillage_post);
}

/*----------------------------------------------------------------------------
 *  Affiche les temps d'exécution
 *----------------------------------------------------------------------------*/

static void
_chrono_total(void)
{
  double  utime;
  double  stime;
  double  time_cpu;
  double  time_tot;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  printf(_("\n\nTime and memory summary\n"
           "-----------------------\n\n")) ;

  ecs_timer_cpu_times(&utime, &stime) ;

  if (utime > 0. || stime > 0.)
    time_cpu = utime + stime;

  else
    time_cpu = ecs_timer_cpu_time() ;

  /* (heure de fin d'execution) - (heure de debut d'execution) */

  if (utime > 0. || stime > 0.) {

    printf("  ") ;
    ecs_print_padded_str(_("User CPU time                       (sec)"),
                         ECS_LNG_AFF_STR) ;
    printf(" : %*.*f\n",
           ECS_LNG_AFF_REE_MANTIS, ECS_LNG_AFF_REE_PRECIS,
           (float)utime);

    printf("  ") ;
    ecs_print_padded_str(_("System CPU time                     (sec)"),
                         ECS_LNG_AFF_STR) ;
    printf(" : %*.*f\n",
           ECS_LNG_AFF_REE_MANTIS, ECS_LNG_AFF_REE_PRECIS,
           (float)stime);
  }

  else if (time_cpu > 0.) {

    printf("  ") ;
    ecs_print_padded_str(_("Total CPU time                      (sec)"),
                         ECS_LNG_AFF_STR) ;
    printf(" : %*.*f\n",
           ECS_LNG_AFF_REE_MANTIS, ECS_LNG_AFF_REE_PRECIS,
           (float)time_cpu);
  }

  /* Durée d'exécution  */

  time_tot = ecs_timer_wtime();

  if (time_tot > 0.) {

    printf("  ") ;
    ecs_print_padded_str(_("Total time                          (sec)"),
                         ECS_LNG_AFF_STR) ;
    printf(" : %*.*f\n",
           ECS_LNG_AFF_REE_MANTIS, ECS_LNG_AFF_REE_PRECIS,
           (float)time_tot);

    if (time_cpu > 0.) {

      printf("  ") ;
      ecs_print_padded_str(_("Total CPU time / Total time              "),
                           ECS_LNG_AFF_STR) ;
      printf(" : %*.*f\n",
             ECS_LNG_AFF_REE_MANTIS, ECS_LNG_AFF_REE_PRECIS,
             (float)(time_cpu/time_tot));

    }

  }
}

/*============================================================================
 *                             Fonction principale
 *============================================================================*/

int
main(int    argc,
     char  *argv[])
{
  ecs_cmd_t         *cmd;

  /* Communication avec le noyau */

  bool               passe_verif;

  ecs_maillage_t    *maillage;

  ecs_post_t        *cas_post;

  ecs_tab_int_t      liste_cel_err;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  ecs_init_gestion_erreur();

#if defined(ENABLE_NLS)

  if (getenv("LANG") != NULL)
     setlocale(LC_ALL,"");
  else
     setlocale(LC_ALL,"C");
  setlocale(LC_NUMERIC,"C");

  bindtextdomain (PACKAGE, LOCALEDIR);
  textdomain(PACKAGE);

#endif

  /* Initialisation comptage et gestion mémoire */

  ecs_mem_usage_init();

  ecs_mem_init(getenv("CS_PREPROCESS_MEM_LOG"));


  /* Lecture de la ligne de commande */

  cmd = ecs_cmd__lit_arg(argc, argv);


  /*==========================================================================*/
  /* Lecture des fichiers d'entrée contenant les maillages                    */
  /*==========================================================================*/

  maillage = _lit_maillage(cmd);

  ecs_maillage__calc_coo_ext(maillage);


  /*==========================================================================*/
  /* Tri des types géométriques pour post-traitement                          */
  /*==========================================================================*/

  ecs_maillage__trie_typ_geo(maillage);
  printf("\n");

  /*========================================================================*/
  /* Passages des couleurs et groupes aux familles                          */
  /* (se fait dès que tous les éléments pouvant porter des familles sont    */
  /*  fusionnés)                                                            */
  /*========================================================================*/

  printf(_("\n\n"
           "Defining families\n"
           "-----------------\n\n"));

  ecs_maillage__cree_famille(maillage);


  /*==========================================================================*/
  /* Préparation du Post-traitement                                           */
  /*==========================================================================*/

  cas_post = _init_post(cmd);


  /*==========================================================================*/
  /* Vérification et correction éventuelle de l'orientation des éléments      */
  /*==========================================================================*/

  liste_cel_err.nbr     = 0;
  liste_cel_err.val     = NULL;

  ecs_maillage__orient_nodal(maillage,
                             &liste_cel_err,
                             cmd->correct_orient);


  /*==========================================================================*/
  /* Écriture de la connectivité nodale sur fichier pour Post-traitement      */
  /*==========================================================================*/

  if (cas_post != NULL) {

    ecs_maillage_post__ecr(_("Fluid Domain"),
                           maillage,
                           ECS_POST_TYPE_VOLUME,
                           cas_post);

    if (liste_cel_err.nbr > 0)
      _post_ele_liste(maillage,
                      liste_cel_err,
                      _("Orientation Error"),
                      ECS_POST_TYPE_ERREUR,
                      cas_post);

  }


  /* On libère les listes des éléments avec problème d'orientation */

  if (liste_cel_err.nbr > 0) {
    ECS_FREE(liste_cel_err.val);
    liste_cel_err.nbr = 0;
  }


  /*==========================================================================*/
  /* Création de la table de connectivité descendante                         */
  /*==========================================================================*/

  /* Passage en connectivité descendante */

  ecs_maillage__connect_descend(maillage);

  printf (_("\nEnd of conversion to descending connectivity"
            "\n--------------------------------------------\n"));

  _aff_taille_maillage(maillage);


  /*========================================================================*/
  /* Vérification du maillage                                               */
  /*========================================================================*/

  if (ecs_maillage__ret_entmail_max(maillage) == ECS_ENTMAIL_CEL) {

    passe_verif = ecs_maillage__verif(maillage, cas_post);

    if (passe_verif == false) {

      ecs_warn();
      printf(_("The mesh has face -> cell connectivity errors.\n"
               "We can go no further."));

    }

  }
  else {

    passe_verif = false;

    ecs_warn();
    printf(_("The mesh does not contain volume elements.\n"
             "We can go no further."));

  }


  /* Tous les post-traitements sont terminés à ce stade */

  ecs_post__detruit_cas(cas_post);


  /*========================================================================*/
  /* Envoi des informations pour le noyau                                   */
  /*========================================================================*/

  if (passe_verif == true)
    ecs_maillage_ncs__ecr(cmd->nom_out, maillage);


  /*==========================================================================*/
  /* Libération en mémoire des structures de maillage                         */
  /*==========================================================================*/

  ecs_maillage__detruit(&maillage);

  ecs_cmd__detruit(cmd);

  _chrono_total();

  /* Bilan mémoire */

  {
    int    nb_div;
    double f_size;
    size_t pr_size;
    static char  unit_prefix[] = {'K', 'M', 'G'};

    printf(_("\nMemory use summary:\n\n"));

    pr_size = ecs_mem_usage_max_pr_size();

    if (pr_size > 0) {
      f_size = pr_size;
      for (nb_div = 0; f_size > 1024. && nb_div < 2; nb_div++)
        f_size /= 1024.;
      printf(_("  Total memory used:                         %15.3f %cb\n"),
             f_size, unit_prefix[nb_div]);
    }

    f_size = ecs_mem_size_max();
    for (nb_div = 0; f_size > 1024. && nb_div < 2; nb_div++)
      f_size /= 1024.;
    printf(_("  Theoretical instrumented dynamic memory:   %15.3f %cb\n"),
           f_size, unit_prefix[nb_div]);

  }


  ecs_mem_end();
  ecs_mem_usage_end();

  /*---------------------------------------------------------*/
  /* Impression d'un message de fin normale du préprocesseur */
  /*---------------------------------------------------------*/

  printf(_("\n\n"
           "  .-----------------------.\n"
           "  |                       |\n"
           "  |  Preprocessor finish  |\n"
           "  |                       |\n"
           "  `-----------------------'\n"
           "\n\n"));


  if (passe_verif == false)
    return EXIT_FAILURE;

  else
    return EXIT_SUCCESS;

}

/*----------------------------------------------------------------------------*/
