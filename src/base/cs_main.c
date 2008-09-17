/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2008 EDF S.A., France
 *
 *     contact: saturne-support@edf.fr
 *
 *     The Code_Saturne Kernel is free software; you can redistribute it
 *     and/or modify it under the terms of the GNU General Public License
 *     as published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version.
 *
 *     The Code_Saturne Kernel is distributed in the hope that it will be
 *     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with the Code_Saturne Kernel; if not, write to the
 *     Free Software Foundation, Inc.,
 *     51 Franklin St, Fifth Floor,
 *     Boston, MA  02110-1301  USA
 *
 *============================================================================*/

/*============================================================================
 * Main program
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <errno.h>
#include <locale.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_config.h>
#include <bft_mem.h>
#include <bft_printf.h>
#include <bft_fp_trap.h>
#include <bft_timer.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <fvm_selector.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_benchmark.h"
#include "cs_couplage.h"
#include "cs_ecs_messages.h"
#include "cs_gui.h"
#include "cs_io.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_mesh_solcom.h"
#include "cs_mesh_quality.h"
#include "cs_mesh_warping.h"
#include "cs_mesh_coherency.h"
#include "cs_multigrid.h"
#include "cs_opts.h"
#include "cs_post.h"
#include "cs_proxy_comm.h"
#include "cs_renumber.h"
#include "cs_sles.h"
#include "cs_suite.h"
#include "cs_syr3_coupling.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * SUBROUTINE CSINIT : sous-programme d'initialisation Fortran listing
 *----------------------------------------------------------------------------*/

extern void CS_PROCF(csinit, CSINIT)
(
 cs_int_t  *ifoenv,    /* Maillage SolCom ou Préprocesseur                    */
 cs_int_t  *iparal,    /* Rang du noyau en cas de parallelisme                */
 cs_int_t  *nparal,    /* Nombre de processus (=1 en sequentiel)              */
 cs_int_t  *ilisr0,    /* Option de sortie du listing principal (rang 0) :    */
                       /*   0 : non redirigé ; 1 : sur fichier 'listing'      */
 cs_int_t  *ilisrp     /* Option de sortie des listings de rang > 0 :         */
                       /*   0 : non redirigé ; 1 : sur fichier 'listing_n*' ; */
                       /*   2 : sur fichier '/dev/null' (suppression)         */
);

/*----------------------------------------------------------------------------
 * SUBROUTINE INITI1 : sous-programme d'initialisation Fortran
 *----------------------------------------------------------------------------*/

extern void CS_PROCF(initi1, INITI1)
(
 cs_int_t  *iverif     /* Activation des tests élémentaires                   */
);

/*----------------------------------------------------------------------------
 * SUBROUTINE CALTRI : sous-programme principal Fortran
 *----------------------------------------------------------------------------*/

extern void CS_PROCF(caltri, CALTRI)
(
 cs_int_t   *iverif,   /* Activation des tests elementaires                   */
 cs_int_t   *nideve,   /* Longueur du tableau d'entiers IDEVEL                */
 cs_int_t   *nrdeve,   /* Longueur du tableau de reels  RDEVEL                */
 cs_int_t   *nituse,   /* Longueur du tableau d'entiers ITUSER                */
 cs_int_t   *nrtuse,   /* Longueur du tableau de reels  RTUSER                */
 cs_int_t   *ifacel,   /* Éléments voisins d'une face interne                 */
 cs_int_t   *ifabor,   /* Élément  voisin  d'une face de bord                 */
 cs_int_t   *ifmfbr,   /* Numéro de famille d'une face de bord                */
 cs_int_t   *ifmcel,   /* Numéro de famille d'une cellule                     */
 cs_int_t   *iprfml,   /* Propriétés d'une famille                            */
 cs_int_t   *ipnfac,   /* Pointeur par sommet dans NODFAC (optionnel)         */
 cs_int_t   *nodfac,   /* Connectivité faces internes/sommets (optionnelle)   */
 cs_int_t   *ipnfbr,   /* Pointeur par sommet dans NODFBR (optionnel)         */
 cs_int_t   *nodfbr,   /* Connectivité faces de bord/sommets (optionnelle)    */
 cs_int_t   *idevel,   /* Pointeur sur le tableau d'entiers IDEVEL            */
 cs_int_t   *ituser,   /* Pointeur sur le tableau d'entiers ITUSER            */
 cs_int_t   *ia,       /* Pointeur sur le tableau d'entiers IA                */
 cs_real_t  *xyzcen,   /* Points associés aux centres des volumes de contrôle */
 cs_real_t  *surfac,   /* Vecteurs surfaces des faces internes                */
 cs_real_t  *surfbo,   /* Vecteurs surfaces des faces de bord                 */
 cs_real_t  *cdgfac,   /* Centres de gravité des faces internes               */
 cs_real_t  *cdgfbr,   /* Centres de gravité des faces de bord                */
 cs_real_t  *xyznod,   /* Coordonnées des sommets (optionnelle)               */
 cs_real_t  *volume,   /* Volumes des cellules                                */
 cs_real_t  *rdevel,   /* Pointeur sur le tableau de reels RDEVEL             */
 cs_real_t  *rtuser,   /* Pointeur sur le tableau de reels RTUSER             */
 cs_real_t  *ra        /* Pointeur sur le tableau de reels RA                 */
);

/*----------------------------------------------------------------------------
 * Fonction utilisateur pour la modification de la géométrie
 *----------------------------------------------------------------------------*/

void CS_PROCF (usmodg, USMODG)
(
 const cs_int_t  *ndim,      /* --> dimension de l'espace                     */
 const cs_int_t  *ncelet,    /* --> nombre de cellules étendu                 */
 const cs_int_t  *ncel,      /* --> nombre de cellules                        */
 const cs_int_t  *nfac,      /* --> nombre de faces internes                  */
 const cs_int_t  *nfabor,    /* --> nombre de faces de bord                   */
 const cs_int_t  *nfml,      /* --> nombre de familles                        */
 const cs_int_t  *nprfml,    /* --> nombre de proprietes des familles         */
 const cs_int_t  *nnod,      /* --> nombre de sommets                         */
 const cs_int_t  *lndfac,    /* --> longueur de nodfac                        */
 const cs_int_t  *lndfbr,    /* --> longueur de nodfbr                        */
 const cs_int_t   ifacel[],  /* --> connectivité faces internes / cellules    */
 const cs_int_t   ifabor[],  /* --> connectivité faces de bord / cellules     */
 const cs_int_t   ifmfbr[],  /* --> liste des familles des faces de bord      */
 const cs_int_t   ifmcel[],  /* --> liste des familles des cellules           */
 const cs_int_t   iprfml[],  /* --> liste des propriétés des familles         */
 const cs_int_t   ipnfac[],  /* --> rang dans nodfac 1er sommet faces int.    */
 const cs_int_t   nodfac[],  /* --> numéro des sommets des faces intérieures  */
 const cs_int_t   ipnfbr[],  /* --> rang dans nodfbr 1er sommet faces bord    */
 const cs_int_t   nodfbr[],  /* --> numéro des sommets des faces de bord      */
       cs_real_t  xyznod[]   /* --> coordonnées des sommets                   */
);

/*----------------------------------------------------------------------------
 * SUBROUTINE MAJGEO : sous-progamme de mise a jour des dimensions du maillage
 *                     dans les commons Fortran
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (majgeo, MAJGEO)
(
 const cs_int_t   *const ncel,    /* --> Number of halo cells             */
 const cs_int_t   *const ncelet,  /* --> Number of halo cells             */
 const cs_int_t   *const nfac,    /* --> Number of internal faces         */
 const cs_int_t   *const nfabor,  /* --> Number of border faces           */
 const cs_int_t   *const nsom,    /* --> Number of vertices               */
 const cs_int_t   *const lndfac,  /* --> Internal face -> vtx array size  */
 const cs_int_t   *const lndfbr,  /* --> Boundary face -> vtx array size  */
 const cs_int_t   *const ncelgb,  /* --> Global number of cells           */
 const cs_int_t   *const nfacgb,  /* --> Global number of internal faces  */
 const cs_int_t   *const nfbrgb,  /* --> Global number of boundary faces  */
 const cs_int_t   *const nsomgb   /* --> Global number of vertices        */
);

/*----------------------------------------------------------------------------
 * SUBROUTINE MEMINI : sous-progamme d'initialisation memoire Fortran
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (memini, MEMINI)
(
 cs_int_t  *iasize,    /* Longueur du tableau d'entiers IA                    */
 cs_int_t  *rasize,    /* Longueur du tableau de reels  RA                    */
 cs_int_t  *nideve,    /* Longueur du tableau d'entiers IDEVEL                */
 cs_int_t  *nrdeve,    /* Longueur du tableau de reels  RDEVEL                */
 cs_int_t  *nituse,    /* Longueur du tableau d'entiers ITUSER                */
 cs_int_t  *nrtuse     /* Longueur du tableau de reels  RTUSER                */
);

/*============================================================================
 * Prototypes de fonctions privées
 *============================================================================*/

/*============================================================================
 * Programme principal
 *============================================================================*/

int main
(
 int    argc,       /* Nombre d'arguments dans la ligne de commandes */
 char  *argv[]      /* Tableau des arguments de la ligne de commandes */
)
{
  double  t1, t2;

  cs_int_t  iasize, rasize;
  cs_int_t  nituse, nrtuse, nideve, nrdeve;

  cs_opts_t  opts;

  int  rang_deb = -1;
  int  _verif = -1;

  cs_int_t  *ia = NULL;
  cs_int_t  *ituser = NULL;
  cs_int_t  *idevel = NULL;

  cs_real_t  *ra = NULL;
  cs_real_t  *rtuser = NULL;
  cs_real_t  *rdevel = NULL;

  /* Première analyse de la ligne de commande pour savoir si l'on a besoin
     de MPI ou non, et initialisation de MPI le cas échéant */

#if defined(_CS_HAVE_MPI)
  rang_deb = cs_opts_mpi_rank(&argc, &argv);
  if (rang_deb > -1)
    cs_base_mpi_init(&argc, &argv, rang_deb);
#endif

  /* initialisation par défaut */

#if defined(_CS_ARCH_Linux)

  if (getenv("LANG") != NULL)
    setlocale(LC_ALL,"");
  else
    setlocale(LC_ALL, "C");
  setlocale(LC_NUMERIC, "C");

#endif

#if defined(ENABLE_NLS)
  bindtextdomain(PACKAGE, LOCALEDIR);
  textdomain(PACKAGE);
#endif

  (void)bft_timer_wtime();

  bft_fp_trap_set();

  /* Initialisation de la gestion mémoire et des signaux */

  cs_base_mem_init();
  cs_base_erreur_init();

  /* interprétation des arguments de la ligne de commande */

  cs_opts_define(argc, argv, &opts);

  /* Ouverture des fichiers 'listing' pour les noeuds de rang > 0 */

  CS_PROCF(csinit, CSINIT)(&(opts.ifoenv),
                           &cs_glob_base_rang,
                           &cs_glob_base_nbr,
                           &(opts.ilisr0),
                           &(opts.ilisrp));
  cs_base_bft_printf_set();

  /* Entête et rappel des options de la ligne de commande */

  cs_opts_logfile_head(argc, argv);

  /* Connexion éventuelle avec le lanceur CFD_Proxy */

  if (opts.proxy_socket != NULL) {
    cs_proxy_comm_initialize(opts.proxy_socket,
                             opts.proxy_key,
                             CS_PROXY_COMM_TYPE_SOCKET);
    BFT_FREE(opts.proxy_socket);
    opts.proxy_key = -1;
  }

  /* Infos système */

  cs_base_info_systeme();

  /* Initialisation des structures globales liées au maillage principal */

  cs_glob_mesh = cs_mesh_create();
  cs_glob_mesh_builder = cs_mesh_builder_create();
  cs_glob_mesh_quantities = cs_mesh_quantities_create();

  /* Initialisation de la lecture des données Préprocesseur */

  if (opts.ifoenv != 0) {

#if defined(FVM_HAVE_MPI)
    cs_glob_pp_io = cs_io_initialize("preprocessor_output",
                                     "Face-based mesh definition, R0",
                                     CS_IO_MODE_READ,
                                     0,
                                     CS_IO_ECHO_OPEN_CLOSE,
                                     cs_glob_base_mpi_comm);
#else
    cs_glob_pp_io = cs_io_initialize("preprocessor_output",
                                     "Face-based mesh definition, R0",
                                     CS_IO_MODE_READ,
                                     CS_IO_ECHO_OPEN_CLOSE,
                                     -1);
#endif

    /* Initialisation des communications avec Syrthes */

    if (cs_syr3_coupling_n_couplings() != 0) {

      cs_int_t coupl_id;
      cs_int_t n_coupl = cs_syr3_coupling_n_couplings();

      for (coupl_id = 0; coupl_id < n_coupl; coupl_id++)
        cs_syr3_coupling_init_comm(cs_syr3_coupling_by_id(coupl_id),
                                   coupl_id + 1,
                                   opts.echo_comm);

    } /* Couplage Syrthes */

  } /* Si ifoenv != 0 */

  /* Allocation de structures internes de l'API F77 pour fichiers suite */

  cs_suite_f77_api_init();

  /* Appel du sous-programme d'initalisation ou de l'aide */

  _verif = opts.iverif;
  if (opts.benchmark > 0 && _verif < 0)
    _verif = 0;

  CS_PROCF(initi1, INITI1)(&_verif);

  if (opts.ifoenv == 0) {

    /* Lecture du fichier au format "SolCom" */

    cs_maillage_solcom_lit(cs_glob_mesh,
                           cs_glob_mesh_quantities);

  }
  else {

    /* Lecture des données issues du Préprocesseur */

    cs_ecs_messages_read_data(cs_glob_mesh,
                              cs_glob_mesh_builder);

  } /* End if ifoenv != 0 */

  /* Initialisation des cas du post-traitement principal */

  cs_post_init_pcp_writer();

  /* Initialisations liées à la construction des halos */

  cs_mesh_init_halo(cs_glob_mesh);

  /* Initialisations liées au parallélisme */

  cs_mesh_init_parall(cs_glob_mesh);

  /* Modification éventuelle de la géométrie */

  CS_PROCF (usmodg, USMODG)(&(cs_glob_mesh->dim),
                            &(cs_glob_mesh->n_cells_with_ghosts),
                            &(cs_glob_mesh->n_cells),
                            &(cs_glob_mesh->n_i_faces),
                            &(cs_glob_mesh->n_b_faces),
                            &(cs_glob_mesh->n_families),
                            &(cs_glob_mesh->n_max_family_items),
                            &(cs_glob_mesh->n_vertices),
                            &(cs_glob_mesh->i_face_vtx_connect_size),
                            &(cs_glob_mesh->b_face_vtx_connect_size),
                            cs_glob_mesh->i_face_cells,
                            cs_glob_mesh->b_face_cells,
                            cs_glob_mesh->b_face_family,
                            cs_glob_mesh->cell_family,
                            cs_glob_mesh->family_item,
                            cs_glob_mesh->i_face_vtx_idx,
                            cs_glob_mesh->i_face_vtx_lst,
                            cs_glob_mesh->b_face_vtx_idx,
                            cs_glob_mesh->b_face_vtx_lst,
                            cs_glob_mesh->vtx_coord);

  /* Découpage des faces "gauche" si nécessaire */

  if (opts.cwf == true) {

    t1 = bft_timer_wtime();
    cs_mesh_warping_cut_faces(cs_glob_mesh, opts.cwf_criterion, opts.cwf_post);
    t2 = bft_timer_wtime();

    bft_printf(_("\n Cutting warped faces (%.3g s)\n"), t2-t1);

  }

  /* Renumérotation en fonction des options du code */

  bft_printf(_("\n Renumbering mesh:\n"));
  bft_printf_flush();
  cs_renumber_mesh(cs_glob_mesh,
                   cs_glob_mesh_quantities);

  /* Initialisation des maillages du post-traitement principal */

  cs_post_init_pcp_maillages();

  /* Mise à jour de certaines dimensions du maillage */

  {
    cs_int_t  n_g_cells, n_g_i_faces, n_g_b_faces, n_g_vertices;

    n_g_cells = cs_glob_mesh->n_g_cells;
    n_g_i_faces = cs_glob_mesh->n_g_i_faces;
    n_g_b_faces = cs_glob_mesh->n_g_b_faces;
    n_g_vertices = cs_glob_mesh->n_g_vertices;

    CS_PROCF (majgeo, MAJGEO)(&(cs_glob_mesh->n_cells),
                              &(cs_glob_mesh->n_cells_with_ghosts),
                              &(cs_glob_mesh->n_i_faces),
                              &(cs_glob_mesh->n_b_faces),
                              &(cs_glob_mesh->n_vertices),
                              &(cs_glob_mesh->i_face_vtx_connect_size),
                              &(cs_glob_mesh->b_face_vtx_connect_size),
                              &n_g_cells,
                              &n_g_i_faces,
                              &n_g_b_faces,
                              &n_g_vertices);
  }

  cs_mesh_print_info(cs_glob_mesh);

  /* Destruction du la structure temporaire servant à la construction du
     maillage principal */

  cs_glob_mesh_builder = cs_mesh_builder_destroy(cs_glob_mesh_builder);

  /* Calcul des grandeurs géométriques associées au maillage */

  bft_printf_flush();

  t1 = bft_timer_wtime();
  cs_mesh_quantities_compute(cs_glob_mesh, cs_glob_mesh_quantities);
  t2 = bft_timer_wtime();

  bft_printf(_("\n Computing geometric quantities (%.3g s)\n"), t2-t1);

  cs_mesh_info(cs_glob_mesh);

  /* Initialisation de la partie selector de la structure maillage */
  cs_mesh_init_selectors();

#if 0
  /* For debugging purposes */
  cs_mesh_dump(cs_glob_mesh);
  cs_mesh_quantities_dump(cs_glob_mesh, cs_glob_mesh_quantities);
#endif

  /* Boucle en temps ou critères de qualité selon options de vérification */

  if (opts.iverif == 0) {
    bft_printf(_("\n Computing quality criteria\n"));
    cs_mesh_quality(cs_glob_mesh, cs_glob_mesh_quantities);
  }

  if (opts.iverif >= 0)
    cs_mesh_coherency_check();

  if (opts.benchmark > 0) {
    int mpi_trace_mode = (opts.benchmark == 2) ? 1 : 0;
    cs_benchmark(mpi_trace_mode);
  }

  if (opts.iverif != 0 && opts.benchmark <= 0) {

    /* Allocation des tableaux de travail */

    CS_PROCF(memini, MEMINI)(&iasize, &rasize,
                             &nideve, &nrdeve, &nituse, &nrtuse);

    bft_printf(_("\n"
                 " --- Main Fortran work arrays:\n"
                 "       LONGIA =   %10d (Number of integers)\n"
                 "       LONGRA =   %10d (Number of reals)\n"
                 "       (%d bytes/integer, %d bytes/real)\n"),
               iasize, rasize,
               sizeof(cs_int_t)/sizeof(char),
               sizeof(cs_real_t)/sizeof(char));

    if (nideve > 0 || nrdeve >0)
      bft_printf(_("\n"
                   " --- Developer Fortran work arrays:\n"
                   "       NIDEVE =   %10d (Number of integer)\n"
                   "       NRDEVE =   %10d (Number of reals)\n"),
                 nideve, nrdeve);

    bft_printf(_("\n"
                 " --- User Fortran work arrays:\n"
                 "       NITUSE =   %10d (Number of integers)\n"
                 "       NRTUSE =   %10d (Number of reals)\n\n"),
               nituse, nrtuse);

    cs_base_mem_init_work(iasize, rasize, &ia, &ra);

    BFT_MALLOC(ituser, nituse, cs_int_t);
    BFT_MALLOC(rtuser, nrtuse, cs_real_t);

    BFT_MALLOC(idevel, nideve, cs_int_t);
    BFT_MALLOC(rdevel, nrdeve, cs_real_t);

    /* Initialisation de la résolution des systèmes linéaires */

    cs_sles_initialize();
    cs_multigrid_initialize();

    /*------------------------------------------------------------------------
     *  appel du sous-programme de gestion de calcul (noyau du code)
     *------------------------------------------------------------------------*/

    CS_PROCF(caltri, CALTRI)(&(opts.iverif),
                             &nideve, &nrdeve, &nituse, &nrtuse,
                             cs_glob_mesh->i_face_cells,
                             cs_glob_mesh->b_face_cells,
                             cs_glob_mesh->b_face_family,
                             cs_glob_mesh->cell_family,
                             cs_glob_mesh->family_item,
                             cs_glob_mesh->i_face_vtx_idx,
                             cs_glob_mesh->i_face_vtx_lst,
                             cs_glob_mesh->b_face_vtx_idx,
                             cs_glob_mesh->b_face_vtx_lst,
                             idevel, ituser, ia,
                             cs_glob_mesh_quantities->cell_cen,
                             cs_glob_mesh_quantities->i_face_normal,
                             cs_glob_mesh_quantities->b_face_normal,
                             cs_glob_mesh_quantities->i_face_cog,
                             cs_glob_mesh_quantities->b_face_cog,
                             cs_glob_mesh->vtx_coord,
                             cs_glob_mesh_quantities->cell_vol,
                             rdevel, rtuser, ra);

    /* Fin de la résolution des systèmes linéaires */

    cs_multigrid_finalize();
    cs_sles_finalize();

    /* les fichiers listing de noeuds > 0 sont fermés dans caltri. */

    /* Libération des tableaux de travail */
    BFT_FREE(ia);
    BFT_FREE(ra);

    BFT_FREE(ituser);
    BFT_FREE(rtuser);

    BFT_FREE(idevel);
    BFT_FREE(rdevel);

  }

  bft_printf(_("\n Destroying structures and ending computation\n"));
  bft_printf_flush();

  /* Libération de structures internes de l'API F77 pour fichiers suite */

  cs_suite_f77_api_finalize();

  /* Libération de la mémoire éventuellement affectée aux couplages */

  cs_syr3_coupling_all_destroy();
#if defined(_CS_HAVE_MPI)
  cs_couplage_detruit_tout();
#endif

  /* Libération de la mémoire associée aux post-traitements */

  cs_post_detruit();

  /* Libération du maillage principal */

  cs_mesh_quantities_destroy(cs_glob_mesh_quantities);
  cs_mesh_destroy(cs_glob_mesh);

  /* Fin de communication éventuelle avec un proxy */

  cs_proxy_comm_finalize();

  /* Temps CPU et finalisation de la gestion mémoire */

  cs_base_bilan_temps();
  cs_base_mem_fin();

  /* retour */

  cs_exit(EXIT_SUCCESS);

  /* jamais appelé normalement, mais pour éviter un warning de compilation */
  return 0;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
