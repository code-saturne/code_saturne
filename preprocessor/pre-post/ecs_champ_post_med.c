/*============================================================================
 *  Definitions des fonctions
 *   associees a la structure `ecs_champ_t' decrivant un champ
 *   et realisant les sorties au format MED
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

#include "cs_config.h"

#if defined(HAVE_MED)

/*============================================================================
 *                                 Visibilite
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


/*----------------------------------------------------------------------------
 *  Fichiers `include' publics  du  paquetage global "MED"
 *----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif

#undef PACKAGE
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef VERSION

#include <med.h>

#undef PACKAGE
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef VERSION

#ifdef __cplusplus
}
#endif


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage global "Post-Traitement"
 *----------------------------------------------------------------------------*/

#include "ecs_post.h"
#include "ecs_post_med.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles des paquetages visibles
 *----------------------------------------------------------------------------*/

#include "ecs_famille.h"
#include "ecs_famille_priv.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_champ.h"


/*----------------------------------------------------------------------------
 *  Fichier  `include' du  paquetage courant associe au fichier courant
 *----------------------------------------------------------------------------*/

#include "ecs_champ_post_med.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' prives   du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_med_priv.h"
#include "ecs_champ_priv.h"


/*============================================================================
 *                              Fonctions privees
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie un pointeur sur un tableau de type `med_int'
 *   dont les valeurs sont converties si necessaire
 *   a partir du tableau `val_ecs' ayant `nbr_val' valeurs de type `ecs_int_t'
 *
 *  On convertit les premières `pas_med' valeurs toutes les `pas_ecs'
 *   valeurs (utile pour ne sortir que les deux premières coordonnées en 2D)
 *
 *  Si le tableau renvoye a ete alloue, la fonction positionne le booleen
 *   `bool_libere' a true
 *----------------------------------------------------------------------------*/

static med_int *
ecs_loc_champ_post_med__cv_int(ecs_int_t   *val_ecs,
                               size_t       nbr_val,
                               size_t       pas_ecs,
                               size_t       pas_med,
                               bool        *bool_libere)
{
  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(val_ecs != NULL);
  assert(nbr_val != 0   );

  if ((pas_med == pas_ecs) && (sizeof(med_int) == sizeof(ecs_int_t))) {

    /* Toutes les valeurs sont à conserver et */
    /* Les entiers de type `med_int' et     */
    /* les entiers de type `ecs_int_t'       */
    /* sont codés sur le meme nombre d'octets */

    /* Aucune conversion n'est nécessaire */

    *bool_libere = false;

    return (med_int *)val_ecs;

  }
  else {

    size_t     iloc, ipas, ipos, ival;
    size_t     nbr_pas;
    med_int  * val_med;

    /* On effectue la conversion de type pour chaque valeur */

    nbr_pas = nbr_val / pas_ecs;

    ECS_MALLOC(val_med, nbr_pas * pas_med, med_int);

    *bool_libere = true;

    ival = 0;

    for (ipas = 0; ipas < nbr_pas; ipas++) {

      ipos = ipas * pas_ecs;

      for (iloc = 0; iloc < pas_med; iloc++)
        val_med[ival++] = (med_int)val_ecs[ipos++];

    }

    return val_med;
  }
}

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie un pointeur sur un tableau de type `med_float'
 *   dont les valeurs sont converties si necessaire
 *   a partir du tableau `val_ecs' ayant `nbr_val' valeurs de type `double'
 *
 *  On convertit les premières `pas_med' valeurs toutes les `pas_ecs'
 *   valeurs (utile pour ne sortir que les deux premières coordonnées en 2D)
 *
 *  Si le tableau renvoye a ete alloue, la fonction positionne le booleen
 *   `bool_libere' a true
 *----------------------------------------------------------------------------*/

static med_float *
ecs_loc_champ_post_med__cv_real(ecs_coord_t  *val_ecs,
                                size_t        nbr_val,
                                size_t        pas_ecs,
                                size_t        pas_med,
                                bool         *bool_libere)
{
  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(val_ecs != NULL);
  assert(nbr_val != 0   );

  if ((pas_med == pas_ecs) && (sizeof(med_float) == sizeof(ecs_coord_t))) {

    /* Toutes les valeurs sont à conserver et */
    /* Les flottants de type `med_float' et   */
    /* les flottants de type `ecs_coord_t'    */
    /* sont codés sur le meme nombre d'octets */

    /* Aucune conversion n'est nécessaire */

    *bool_libere = false;

    return (med_float *)val_ecs;

  }
  else {

    size_t       iloc, ipas, ipos, ival;
    size_t       nbr_pas;
    med_float  * val_med;

    /* On effectue la conversion de type pour chaque valeur */

    nbr_pas = nbr_val / pas_ecs;

    ECS_MALLOC(val_med, nbr_pas * pas_med, med_float);

    *bool_libere = true;

    ival = 0;

    for (ipas = 0; ipas < nbr_pas; ipas++) {

      ipos = ipas * pas_ecs;

      for (iloc = 0; iloc < pas_med; iloc++)
        val_med[ival++] = (med_float)val_ecs[ipos++];

    }

    return val_med;
  }
}

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie un pointeur sur une sous-structure associée
 *  à un maillage pour un cas de sortie MED.
 *---------------------------------------------------------------------------*/

static ecs_med_maillage_t  *
ecs_loc_champ_post_med__maillage(const ecs_med_t  *cas_med,
                                 const char       *nom_maillage)
{
  ecs_int_t  ind;
  ecs_med_maillage_t  *maillage_med = NULL;

  /* Recherche du maillage */
  /*-----------------------*/

  for (ind = 0; ind < cas_med->nbr_maillages; ind++) {
    maillage_med = cas_med->tab_maillages[ind];
    if (strcmp(nom_maillage, maillage_med->nom_maillage) == 0)
      break;
  }

  if (ind >= cas_med->nbr_maillages)
    ecs_error(__FILE__, __LINE__, 0,
              _("MED: no mesh named \"%s\".\n"
                "is associated with file: \"%s\"\n"),
              nom_maillage, cas_med->nom_fic);

  return  maillage_med;
}

/*----------------------------------------------------------------------------
 *  Fonction écrivant une famille
 *
 *  Les numéros des descripteurs de type "groupe" seront perdus.
 *----------------------------------------------------------------------------*/

static void
ecs_loc_champ_post_med__ecr_fam(const char           *prefixe_nom_fam,
                                char                 *nom_maillage_med,
                                const med_int         num_fam_med,
                                const ecs_famille_t  *ptr_fam,
                                ecs_med_t            *cas_med)
{
  size_t   ipropr;

  ecs_tab_char_t  tab_nom_descr;
  ecs_tab_int_t   tab_ide_descr;

  /* Déclarations des variables pour MED */
  /*-------------------------------------*/

  char        nom_fam_med[MED_TAILLE_NOM + 1];
  char        str_num_fam_med[MED_TAILLE_NOM + 1];
  char      * grp_nom_med;
  char      * ptr_grp_nom_med;

  size_t      ind;
  med_int     nbr_grp_med;

  med_err     ret_med = 0;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


  /* Récupération des propriétés de la famille */

  tab_nom_descr = ecs_famille__ret_nom(ptr_fam);

  nbr_grp_med = (med_int)tab_nom_descr.nbr;


  /* Construction du numéro et du nom de la famille */

  sprintf(str_num_fam_med, "%d", (int)num_fam_med);

  assert (MED_TAILLE_NOM > 3);

  strncpy(nom_fam_med, prefixe_nom_fam, MED_TAILLE_NOM - 3);
  strncat(nom_fam_med, str_num_fam_med,
          MED_TAILLE_NOM - strlen(prefixe_nom_fam));
  nom_fam_med[MED_TAILLE_NOM] = '\0';


  /* Création des groupes MED */

  if (nbr_grp_med > 0)
    ECS_MALLOC(grp_nom_med, MED_TAILLE_LNOM  * nbr_grp_med + 1, char);
  else
    grp_nom_med = NULL;

  /* Affectation des groupes MED */

  for (ipropr = 0, nbr_grp_med = 0;
       ipropr < tab_nom_descr.nbr;
       ipropr++) {

    ptr_grp_nom_med = grp_nom_med + (MED_TAILLE_LNOM*nbr_grp_med);
    ind = 0;
    while (   ind < MED_TAILLE_LNOM
           && tab_nom_descr.val[ipropr][ind] != '\0') {
      ptr_grp_nom_med[ind] = tab_nom_descr.val[ipropr][ind];
      ind++;
    }
    while (ind < MED_TAILLE_LNOM)
      ptr_grp_nom_med[ind++] = ' ';

    nbr_grp_med++;

  }

  ECS_FREE(tab_ide_descr.val);
  ECS_FREE(tab_nom_descr.val);

  /* Appel à la fonction d'écriture MED */

  ret_med = MEDfamCr(cas_med->fid,
                     nom_maillage_med,
                     nom_fam_med,
                     num_fam_med,
                     NULL,
                     NULL,
                     NULL,
                     0,
                     grp_nom_med,
                     nbr_grp_med);

  if (ret_med != 0)
    ecs_error(__FILE__, __LINE__, 0,
              _("MED: error writing file \"%s\".\n"
                "Name   of family to write: \"%s\"\n"
                "Number of family to write: \"%d\""),
              cas_med->nom_fic, nom_fam_med, (int)num_fam_med);

  ECS_FREE(grp_nom_med);
}

/*============================================================================
 *                             Fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction écrivant les familles
 *----------------------------------------------------------------------------*/

void
ecs_champ_post_med__ecr_famille(const char           *nom_maillage,
                                const ecs_famille_t  *famille_elt,
                                const ecs_famille_t  *famille_inf,
                                ecs_med_t            *cas_med)
{
  ecs_int_t        ifam_ent;
  size_t           ind;

  ecs_med_maillage_t  *maillage_med = NULL;

  const ecs_famille_t  *ptr_fam;

  /* Declarations des variables pour MED */
  /*-------------------------------------*/

  char        nom_fam_med[MED_TAILLE_NOM + 1];
  char        str_num_fam_med[MED_TAILLE_NOM + 1];

  med_int     num_fam_med;

  med_err     ret_med = 0;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Recherche du maillage */
  /*-----------------------*/

  for (ind = 0; ind < (size_t)(cas_med->nbr_maillages); ind++) {
    maillage_med = cas_med->tab_maillages[ind];
    if (strcmp(nom_maillage, maillage_med->nom_maillage) == 0)
      break;
  }

  if (ind >= (size_t)(cas_med->nbr_maillages))
    ecs_error(__FILE__, __LINE__, 0,
              _("MED: no mesh named \"%s\".\n"
                "is associated with file: \"%s\"\n"),
              nom_maillage, cas_med->nom_fic);


  /* Création de la famille 0 */
  /*--------------------------*/

  num_fam_med = 0;

  sprintf(str_num_fam_med, "%d", num_fam_med);

  assert(strlen("FAMILLE_") + strlen(str_num_fam_med) <= MED_TAILLE_NOM);

  strcpy(nom_fam_med, "FAMILLE_");
  strcat(nom_fam_med, str_num_fam_med);

  ret_med = MEDfamCr(cas_med->fid,
                     maillage_med->nom_maillage_med,
                     nom_fam_med,
                     num_fam_med,
                     NULL,
                     NULL,
                     NULL,
                     0,
                     NULL,
                     0);

  if (ret_med != 0)
    ecs_error(__FILE__, __LINE__, 0,
              _("MED: error writing file \"%s\".\n"
                "Name   of family to write: \"%s\"\n"
                "Number of family to write: \"%d\""),
              cas_med->nom_fic, nom_fam_med, (int)num_fam_med);

  /* On affecte la famille 0 à tous les sommets */
  /*--------------------------------------------*/

  /* Création des familles MED pour les elements */
  /*---------------------------------------------*/

  for (ifam_ent = 0; ifam_ent < 2; ifam_ent ++) {

    if (ifam_ent == 0)
      ptr_fam = famille_elt;
    else
      ptr_fam = famille_inf;

    for (; ptr_fam != NULL; ptr_fam  = ptr_fam->l_famille_sui) {

      num_fam_med = - (ptr_fam->num);
      ecs_loc_champ_post_med__ecr_fam("FAMILLE_ELEMENT_",
                                      maillage_med->nom_maillage_med,
                                      num_fam_med,
                                      ptr_fam,
                                      cas_med);

    }
  } /* Fin : boucle sur les tetes de liste chaînée des familles */
}

/*----------------------------------------------------------------------------
 *  Fonction imprimant le contenu des tables asociees aux sommets
 *----------------------------------------------------------------------------*/

void
ecs_champ_post_med__ecr_som(const char         *nom_maillage,
                            size_t              n_vertices,
                            ecs_coord_t         vertex_coords[],
                            const ecs_med_t    *cas_med)
{
  bool          bool_libere_coo_noe;

  ecs_int_t     ind;
  ecs_int_t     isom;
  ecs_int_t     nbr_som;

  ecs_med_maillage_t  *maillage_med;

  /* Declarations des variables pour MED */
  /*-------------------------------------*/

  char      *nom_coo_med;
  char      *uni_coo_med;


  med_err    ret_med = 0;

  med_int    nbr_noe_med;

  med_int   *fam_noe_med;

  med_float *coo_noe_med;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Recherche du maillage med */
  /*---------------------------*/

  maillage_med = ecs_loc_champ_post_med__maillage(cas_med,
                                                  nom_maillage);

  /* Nombre de noeuds */
  /*------------------*/

  nbr_som     = n_vertices;
  nbr_noe_med = (med_int)nbr_som;

  /* Noms et unites des coordonnees */
  /*--------------------------------*/

  ECS_MALLOC(nom_coo_med, 3 * MED_TAILLE_PNOM + 1, char);
  ECS_MALLOC(uni_coo_med, 3 * MED_TAILLE_PNOM + 1, char);
  for (ind = 0; ind < (ecs_int_t)(3 * MED_TAILLE_PNOM); ind++)
    nom_coo_med[ind] = ' ', uni_coo_med[ind] = ' ';
  nom_coo_med[MED_TAILLE_PNOM * 3] = '\0';
  uni_coo_med[MED_TAILLE_PNOM * 3] = '\0';

  nom_coo_med[0] = 'x';
  nom_coo_med[MED_TAILLE_PNOM    ] = 'y';
  nom_coo_med[MED_TAILLE_PNOM * 2] = 'z';

  /* Familles MED */
  /*--------------*/

  /* On attribue a tous les sommets la famille `0' */

  ECS_MALLOC(fam_noe_med, nbr_som, med_int);

  for (isom = 0; isom < nbr_som; isom++)
    fam_noe_med[isom] = 0;

  /* Coordonnées des noeuds */
  /*------------------------*/

  coo_noe_med = ecs_loc_champ_post_med__cv_real(vertex_coords,
                                                n_vertices,
                                                3,
                                                3,
                                                &bool_libere_coo_noe);

  ret_med = MEDnoeudsEcr(cas_med->fid,
                         maillage_med->nom_maillage_med,
                         (med_int)3,
                         coo_noe_med,
                         MED_FULL_INTERLACE,
                         MED_CART,
                         nom_coo_med,
                         uni_coo_med,
                         NULL,
                         MED_FAUX,
                         NULL,
                         MED_FAUX,
                         fam_noe_med,
                         nbr_noe_med);

  if (ret_med != 0)
    ecs_error(__FILE__, __LINE__, 0,
              _("MED: error writing file \"%s\".\n"
                "Error writing coordinates."),
              cas_med->nom_fic);

  ECS_FREE(nom_coo_med);
  ECS_FREE(uni_coo_med);

  ECS_FREE(fam_noe_med);

  if (bool_libere_coo_noe == true)
    ECS_FREE(coo_noe_med);
}

/*----------------------------------------------------------------------------
 *  Fonction qui écrit les connectivités des éléments
 *   selon leur type géometrique
 *
 *  Les éléments doivent avoir ete triés suivant leur type géometrique
 *----------------------------------------------------------------------------*/

void
ecs_champ_post_med__ecr_elt(const char           *nom_maillage,
                            ecs_champ_t          *champ_def,
                            const int             elt_fam[],
                            const ecs_tab_int_t  *tab_elt_typ_geo,
                            const ecs_med_t      *cas_med)
{
  size_t      cpt_elt;
  int         elt_typ_ref;
  size_t      ind;
  size_t      ielt;
  size_t      ifac;
  size_t      isom;
  size_t      ival;
  int         marqueur_fin;
  size_t      nbr_fac_loc;
  size_t      nbr_elt;
  size_t      nbr_elt_typ_geo;   /* Nb. elements de meme type geometrique */
  size_t      nbr_som_elt;
  size_t      nbr_som_loc;
  size_t      pos_elt;
  size_t      renum_som;

  ecs_size_t  *def_pos_tab;
  ecs_int_t   *def_val_tab;

  ecs_med_maillage_t  * maillage_med;

  /* Declarations des variables pour MED */
  /*-------------------------------------*/

  med_geometrie_element  typ_geo_med;

  med_err            ret_med = 0;

  med_int            ityp_med;
  med_int            mdim_med;
  med_int            nbr_ele_med;
  med_int            nbr_som_med;

  med_int          * index_med    = NULL;
  med_int          * index_f_med  = NULL;
  med_int          * connect_med  = NULL;
  med_int          * fam_ele_med  = NULL;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(champ_def != NULL);

  nbr_elt = champ_def->nbr;

  /* Recherche du maillage med */
  /*---------------------------*/

  maillage_med = ecs_loc_champ_post_med__maillage(cas_med,
                                                  nom_maillage);

  /* Familles MED des elements */
  /*---------------------------*/

  ECS_MALLOC(fam_ele_med, nbr_elt, med_int);

  if (elt_fam != NULL) {
    for (ielt = 0; ielt < nbr_elt; ielt++)
      fam_ele_med[ielt] = -elt_fam[ielt];
  }
  else {
    for (ielt = 0; ielt < nbr_elt; ielt++)
      fam_ele_med[ielt] = 0;
  }

  /* Connectivite des elements */
  /*---------------------------*/

  mdim_med = maillage_med->dim_entite;

  ecs_champ__regle_en_pos(champ_def);

  def_pos_tab = champ_def->pos;
  def_val_tab = champ_def->val;

  /* Boucle sur les éléments ayant le même type géométrique */
  /*--------------------------------------------------------*/

#define ECS_FCT_TYP_ECS(ityp_med) ecs_fic_med_init_elt_liste_c[ityp_med].ecs_typ

  cpt_elt = 0;

  elt_typ_ref = -1;

  while (cpt_elt < nbr_elt) {

    /* Recherche du prochain type d'élément utilisé */

    elt_typ_ref += 1;

    while (tab_elt_typ_geo->val[elt_typ_ref] == 0)
      elt_typ_ref++;

    /* détermination du type géométrique MED correspondant */

    ityp_med = 0;
    while (ityp_med < ECS_MED_NBR_TYP_ELT
           && ((ecs_int_t)ecs_fic_med_init_elt_liste_c[ityp_med].ecs_typ
               != elt_typ_ref))
      ityp_med++;

    if (ityp_med == ECS_MED_NBR_TYP_ELT)
      ecs_error(__FILE__, __LINE__, 0,
                _("MED: error writing file \"%s\".\n"
                  "The element geometric type has no MED equivalent.\n"
                  "Element geometric type: \"%d\""),
                cas_med->nom_fic, (int)elt_typ_ref);

    typ_geo_med = ecs_fic_med_init_elt_liste_c[ityp_med].med_typ;

    /* On compte le nombre d'éléments ayant le même type géométrique */

    nbr_elt_typ_geo = tab_elt_typ_geo->val[elt_typ_ref];

    pos_elt = def_pos_tab[cpt_elt] - 1;

    nbr_som_elt
      = ecs_fic_elt_typ_liste_c[ECS_FCT_TYP_ECS(ityp_med)].nbr_som;

    /* Cas de éléments "classiques" (non polygonaux/polyédriques) */
    /*------------------------------------------------------------*/

    if (typ_geo_med != MED_POLYGONE && typ_geo_med != MED_POLYEDRE) {

      /* Prise en compte des définitions des connectivités MED */
      /*-------------------------------------------------------*/

      ECS_MALLOC(connect_med, nbr_elt_typ_geo * nbr_som_elt, med_int);

      for (ielt = 0; ielt < nbr_elt_typ_geo; ielt++) {

        for (isom = 0; isom < nbr_som_elt; isom++) {

          renum_som = ecs_fic_med_init_elt_liste_c[ityp_med].num_som[isom] - 1;

          connect_med[ielt * nbr_som_elt + isom]
            = def_val_tab[pos_elt + (ielt * nbr_som_elt) + renum_som];

        }
      }

      /* Convention sur la taille des mailles */

      nbr_ele_med = (med_int)nbr_elt_typ_geo;

      ret_med = MEDconnEcr(cas_med->fid,
                           maillage_med->nom_maillage_med,
                           mdim_med,
                           connect_med,
                           MED_FULL_INTERLACE,
                           nbr_ele_med,
                           MED_MAILLE,
                           typ_geo_med,
                           MED_NOD);

      if (ret_med != 0)
        ecs_error(__FILE__, __LINE__, 0,
                  _("MED: error writing file \"%s\".\n"
                    "Error writing connectivity."),
                  cas_med->nom_fic);

      ECS_FREE(connect_med);
    }

    /* Cas de éléments polygonaux */
    /*----------------------------*/

    else if (typ_geo_med == MED_POLYGONE) {

      nbr_som_med =   def_pos_tab[cpt_elt + nbr_elt_typ_geo]
                    - def_pos_tab[cpt_elt];

      /* Recopie (avec translation d'index) des définitions */

      ECS_MALLOC(index_med, nbr_elt_typ_geo + 1, med_int);
      ECS_MALLOC(connect_med, nbr_som_med, med_int);

      for (ind = 0; ind < nbr_elt_typ_geo + 1; ind++)
        index_med[ind] = def_pos_tab[cpt_elt + ind] - pos_elt;

      for (ind = 0; ind < (size_t)nbr_som_med; ind++)
        connect_med[ind] = def_val_tab[pos_elt + ind];

      ret_med = MEDpolygoneConnEcr(cas_med->fid,
                                   maillage_med->nom_maillage_med,
                                   index_med,
                                   (med_int)(nbr_elt_typ_geo + 1),
                                   connect_med,
                                   MED_MAILLE,
                                   MED_NOD);

      if (ret_med != 0)
        ecs_error(__FILE__, __LINE__, 0,
                  _("MED: error writing file \"%s\".\n"
                    "Error writing polygons connectivity."),
                  cas_med->nom_fic);

      ECS_FREE(index_med);
      ECS_FREE(connect_med);
    }

    /* Cas de éléments polyèdriques */
    /*------------------------------*/

    else if (typ_geo_med == MED_POLYEDRE) {

      /* Convention : définition nodale cellule->sommets avec numéros de
         premiers sommets répétés en fin de liste pour marquer la fin
         de chaque face */

      /* Index éléments -> faces */

      ECS_MALLOC(index_med, nbr_elt_typ_geo + 1, med_int);

      index_med[0] = 1;

      for (ielt = 0; ielt < nbr_elt_typ_geo; ielt++) {

        marqueur_fin = -1;
        nbr_fac_loc = 0;

        for (ival = def_pos_tab[cpt_elt + ielt    ] - 1;
             ival < def_pos_tab[cpt_elt + ielt + 1] - 1;
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

        index_med[ielt + 1] = index_med[ielt] + nbr_fac_loc;
      }

      /* Index faces -> sommets */

      ECS_MALLOC(index_f_med, index_med[nbr_elt_typ_geo], med_int);

      ifac = 0;
      marqueur_fin = -1;
      nbr_som_loc = 0;

      index_f_med[0] = 1;

      for (ival = def_pos_tab[cpt_elt                  ] - 1;
           ival < def_pos_tab[cpt_elt + nbr_elt_typ_geo] - 1;
           ival++) {

        if (def_val_tab[ival] != marqueur_fin) {
          nbr_som_loc += 1;
          if (marqueur_fin == -1)
            marqueur_fin = def_val_tab[ival];
        }
        else {
          index_f_med[ifac + 1] = index_f_med[ifac] + nbr_som_loc;
          ifac += 1;
          marqueur_fin = -1;
          nbr_som_loc = 0;
        }
      }

      assert(ifac == (size_t)(index_med[nbr_elt_typ_geo] - 1));

      /* Connectivité */

      nbr_som_med =  index_f_med[index_med[nbr_elt_typ_geo] - 1] - 1;

      ECS_MALLOC(connect_med, nbr_som_med, med_int);

      isom = 0;
      marqueur_fin = -1;

      for (ival = def_pos_tab[cpt_elt                  ] - 1;
           ival < def_pos_tab[cpt_elt + nbr_elt_typ_geo] - 1;
           ival++) {

        if (def_val_tab[ival] != marqueur_fin) {
          connect_med[isom++] = def_val_tab[ival];
          if (marqueur_fin == -1)
            marqueur_fin = def_val_tab[ival];
        }
        else
          marqueur_fin = -1;

      }

      assert(isom == (size_t)nbr_som_med);

      ret_med = MEDpolyedreConnEcr(cas_med->fid,
                                   maillage_med->nom_maillage_med,
                                   index_med,
                                   (med_int)(nbr_elt_typ_geo + 1),
                                   index_f_med,
                                   index_med[nbr_elt_typ_geo],
                                   connect_med,
                                   MED_NOD);

      if (ret_med != 0)
        ecs_error(__FILE__, __LINE__, 0,
                  _("MED: error writing file \"%s\".\n"
                    "Error writing polyhedra connectivity."),
                  cas_med->nom_fic);


      ECS_FREE(index_med);
      ECS_FREE(index_f_med);
      ECS_FREE(connect_med);
    }

    /* Familles MED des éléments */
    /*---------------------------*/

    if (typ_geo_med != MED_POLYGONE && typ_geo_med != MED_POLYEDRE) {

      ret_med = MEDfamEcr(cas_med->fid,
                          maillage_med->nom_maillage_med,
                          fam_ele_med + cpt_elt,
                          nbr_elt_typ_geo,
                          MED_MAILLE,
                          typ_geo_med);

      if (ret_med != 0)
        ecs_error(__FILE__, __LINE__, 0,
                  _("MED: error writing file \"%s\".\n"
                    "Error writing families."),
                  cas_med->nom_fic);
    }

    cpt_elt += nbr_elt_typ_geo;
  }

#undef ECS_FCT_TYP_ECS

  ECS_FREE(fam_ele_med);

  ecs_champ__libere_pos_tab(champ_def, def_pos_tab);
}

/*----------------------------------------------------------------------------
 *  Fonction qui ajoute à une structure maillage_med les informations
 *   sur le nombre d'éléments de chaque type d'un maillage
 *----------------------------------------------------------------------------*/

void
ecs_champ_post_med__cpt_elt_typ(const ecs_tab_int_t  *tab_elt_typ_geo,
                                const char           *nom_maillage,
                                ecs_med_t            *cas_med)
{
  size_t  cpt_typ_med;
  size_t  ityp;
  int     ityp_med;

  ecs_med_maillage_t  * maillage_med;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Recherche du maillage med */
  /*---------------------------*/

  maillage_med = ecs_loc_champ_post_med__maillage(cas_med,
                                                  nom_maillage);

  if (maillage_med == NULL)
    return;

  /* Comptage et allocation */

  cpt_typ_med = 0;

  for (ityp = 0; ityp < tab_elt_typ_geo->nbr; ityp++) {
    if (tab_elt_typ_geo->val[ityp] > 0)
      cpt_typ_med += 1;
  }

  maillage_med->nbr_typ_ele = cpt_typ_med;
  ECS_MALLOC(maillage_med->nbr_ele_typ, cpt_typ_med, ecs_int_t);
  ECS_MALLOC(maillage_med->med_typ,     cpt_typ_med, med_geometrie_element);

  /* Mise à jour de la structure maillage_med */

  cpt_typ_med = 0;

  for (ityp = 0; ityp < tab_elt_typ_geo->nbr; ityp++) {

    if (tab_elt_typ_geo->val[ityp] > 0) {

      /* détermination du type géométrique MED correspondant */

      ityp_med = 0;
      while (   ityp_med < ECS_MED_NBR_TYP_ELT
             && (   (int)ecs_fic_med_init_elt_liste_c[ityp_med].ecs_typ
                 != (int)ityp))
        ityp_med++;

      if (ityp_med == ECS_MED_NBR_TYP_ELT)
        ecs_error(__FILE__, __LINE__, 0,
                  _("Geometric type \"%d\" has no MED equivalent."),
                  (int)ityp);

      maillage_med->nbr_ele_typ[cpt_typ_med] = tab_elt_typ_geo->val[ityp];
      maillage_med->med_typ[cpt_typ_med]
        = ecs_fic_med_init_elt_liste_c[ityp_med].med_typ;

      cpt_typ_med += 1;

    }
  }
}

/*----------------------------------------------------------------------------
 *  Fonction ecrivant un champ au format MED
 *----------------------------------------------------------------------------*/

void
ecs_champ_post_med__ecr_val(const ecs_tab_int_t  *tab_val,
                            const char           *nom_maillage,
                            const char           *nom_champ,
                            const ecs_med_t      *cas_med)
{
  bool        bool_libere_val;

  ecs_int_t   cpt_elt;
  ecs_int_t   elt_typ_ref;
  ecs_int_t   nbr_elt;
  ecs_int_t   nbr_elt_typ_geo; /* Nombre d'elements de meme type geometrique */

  ecs_med_maillage_t  * maillage_med;

  /* Declarations des variables pour MED */
  /*-------------------------------------*/

  char   *nom_champ_med;
  char    nom_unite_dt_med[MED_TAILLE_NOM + 1] = "";
  char    profil_med_nopfl[] = MED_NOPFL;
  char    locname_med_nogauss[] = MED_NOGAUSS;

  med_geometrie_element  typ_geo_med;

  med_int   ityp_med;
  med_int  *val_med;

  med_err   ret_med = 0;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(tab_val != NULL);

  /* Recherche du maillage med */
  /*---------------------------*/

  maillage_med = ecs_loc_champ_post_med__maillage(cas_med,
                                                  nom_maillage);

  if (maillage_med == NULL)
    return;

  /* Nom du champ */
  /*--------------*/

  ECS_MALLOC(nom_champ_med, strlen(nom_champ) + 1, char);
  strcpy(nom_champ_med, nom_champ);

  /* Valeurs du champ suivant le type des elements */
  /*-----------------------------------------------*/

  nbr_elt = tab_val->nbr;
  cpt_elt = 0;

  elt_typ_ref = -1;

  for (ityp_med = 0; ityp_med < maillage_med->nbr_typ_ele; ityp_med++) {

    typ_geo_med     = maillage_med->med_typ[ityp_med];
    nbr_elt_typ_geo = maillage_med->nbr_ele_typ[ityp_med];

    /* En cas d'incohérence, on sort de la boucle, on gèrera l'erreur ensuite */

    if (cpt_elt + nbr_elt_typ_geo > nbr_elt)
      break;

    if (typ_geo_med != MED_POLYGONE && typ_geo_med != MED_POLYEDRE) {

      /* On écrit les valeurs correspondant à ce type géométrique */

      val_med = ecs_loc_champ_post_med__cv_int(tab_val->val + cpt_elt,
                                               1,
                                               1,
                                               nbr_elt_typ_geo,
                                               &bool_libere_val);

      ret_med = MEDchampEcr(cas_med->fid,
                            maillage_med->nom_maillage_med,
                            nom_champ_med,
                            (unsigned char *)(val_med),
                            MED_FULL_INTERLACE,
                            nbr_elt_typ_geo,
                            locname_med_nogauss,
                            MED_ALL,
                            profil_med_nopfl,
                            MED_NO_PFLMOD,
                            MED_MAILLE,
                            typ_geo_med,
                            MED_NOPDT,
                            nom_unite_dt_med,
                            0.0,
                            MED_NONOR);

      if (ret_med != 0)
        ecs_error(__FILE__, __LINE__, 0,
                  _("MED: error writing file \"%s\".\n"
                    "Error writing field \"%s\"."),
                  cas_med->nom_fic, nom_champ_med);

      if (bool_libere_val == true)
        ECS_FREE(val_med);

    }

    cpt_elt += nbr_elt_typ_geo;
  }

  if (cpt_elt != nbr_elt || ityp_med < maillage_med->nbr_typ_ele)
    ecs_error(__FILE__, __LINE__, 0,
              _("MED: error writing field \"%s\".\n"
                "Incompatibility between the number of elements to write (%d)\n"
                "and the dimensions of mesh \"%s\"."),
              nom_champ_med, (ecs_int_t)nbr_elt, nom_maillage);

  ECS_FREE(nom_champ_med);
}

#endif /* HAVE_MED */

/*----------------------------------------------------------------------------*/

