/*============================================================================
 *  Définition de la fonction
 *   de lecture d'un fichier de maillage au format MED
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

#if defined(HAVE_MED)


/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie standard C
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/

#include "ecs_def.h"
#include "ecs_elt_typ_liste.h"
#include "ecs_mem.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles des paquetages visibles
 *----------------------------------------------------------------------------*/

#include "ecs_descr.h"
#include "ecs_descr_chaine.h"
#include "ecs_famille.h"
#include "ecs_famille_chaine.h"
#include "ecs_maillage.h"
#include "ecs_maillage_priv.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_maillage_pre.h"


/*----------------------------------------------------------------------------
 *  Fichier  `include' du  paquetage courant associe au fichier courant
 *----------------------------------------------------------------------------*/

#include "ecs_pre_med.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' prives   du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_med_priv.h"


/*============================================================================
 *                              Fonctions privees
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Création d'une structure `ecs_med_t' et ouverture d'un fichier MED
 *  en lecture
 *----------------------------------------------------------------------------*/

static ecs_med_t *
ecs_pre_med__cree(const char  *nom_fichier)
{
  ecs_med_t  * fic;

  med_err retval = 0;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  ECS_MALLOC(fic, 1, ecs_med_t);

  ECS_MALLOC(fic->nom_fic, strlen(nom_fichier) + 1, char);
  strcpy(fic->nom_fic, nom_fichier);

  fic->fid = MEDfileOpen(fic->nom_fic, MED_ACC_RDONLY);

  if (fic->fid < 0)
    ecs_error(__FILE__, __LINE__, 0,
              _("MED: error opening file \"%s\"."), nom_fichier);

  retval = MEDfileNumVersionRd(fic->fid,
                               &(fic->version[0]),
                               &(fic->version[1]),
                               &(fic->version[2]));

  if (retval < 0)
    ecs_error(__FILE__, __LINE__, 0,
              _("MED: error checking file version \"%s\"."), nom_fichier);

  fic->nbr_maillages = 0;
  fic->tab_maillages = NULL;

  return fic;
}

/*----------------------------------------------------------------------------
 *  Fermeture d'un fichier MED en lecture et destruction de la structure
 *  `ecs_med_t' associée
 *----------------------------------------------------------------------------*/

static ecs_med_t *
ecs_pre_med__detruit(ecs_med_t  * fic)
{
  med_err retval = 0;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(fic != NULL);

  retval = MEDfileClose(fic->fid);

  if (retval != 0)
    ecs_error(__FILE__, __LINE__, 0,
              _("MED: error closing file \"%s\"."),
              fic->nom_fic);

  ECS_FREE(fic->nom_fic);
  ECS_FREE(fic);

  return fic;
}

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie un pointeur sur un tableau de type `ecs_int_t'
 *   dont les valeurs sont converties si nécessaire
 *   à partir du tableau de connectivité `val_med' ayant `nbr_val' valeurs
 *   de type `med_int', avec `pas_med' valeurs MED et `pas_ecs' valeurs
 *   locales par élement.
 *
 *  Cette fonction supprime donc les références inutiles (noeuds
 *   quadratiques ou connectivités supplémentaires éventuelles)
 *
 *  Le tableau en entrée `val_med' est libéré si le tableau renvoyé est alloué
 *----------------------------------------------------------------------------*/

static ecs_int_t *
ecs_loc_convert_connect_med_ecs(med_int          *val_med,
                                const ecs_int_t   nbr_val,
                                const ecs_int_t   pas_med,
                                const ecs_int_t   pas_ecs)
{
  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(val_med != NULL);
  assert(nbr_val != 0);
  assert(pas_med != 0);
  assert(pas_ecs != 0);

  if ((pas_med == pas_ecs) && (sizeof(med_int) == sizeof(ecs_int_t))) {

    /* Le nombre de références par élément    */
    /* MED est le même que celui par élément  */
    /* local (élément linéaire, pas de        */
    /* connectivité supplémentaire)           */

    /* Les entiers de type `med_int' et       */
    /* les entiers de type `ecs_int_t'        */
    /* sont codés sur le meme nombre d'octets */

    /* Aucune conversion n'est nécessaire     */

    return (ecs_int_t *)val_med;

  }
  else {

    ecs_int_t    ielt;
    ecs_int_t    iloc;
    ecs_int_t    ipos_ecs;
    ecs_int_t    ipos_med;
    ecs_int_t    nbr_val_ecs;
    ecs_int_t  * val_ecs;

    /* Les entiers de type `med_int' et `ecs_int_t'  */
    /* ne sont pas codés sur le même nombre d'octets */

    /* On effectue la conversion de type pour chaque valeur */

    nbr_val_ecs = (nbr_val / pas_med) * pas_ecs;

    ECS_MALLOC(val_ecs, nbr_val_ecs, ecs_int_t);

    ipos_ecs = 0;

    for (ielt = 0; ielt < (nbr_val / pas_med); ielt++) {

      ipos_med = ielt * pas_med;

      for (iloc = 0; iloc < pas_ecs; iloc++)
        val_ecs[ipos_ecs++] = (ecs_int_t)val_med[ipos_med++];

    }

    ECS_FREE(val_med);

    return val_ecs;
  }
}

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie un pointeur sur un tableau de type `double'
 *   dont les valeurs sont converties si nécessaire
 *   à partir du tableau `val_med' ayant `nbr_val' valeurs de type `med_float'
 *
 *  Le tableau en entrée `val_med' est libere si le tableau renvoyé est alloué
 *----------------------------------------------------------------------------*/

static ecs_coord_t *
ecs_loc_convert_real_med_ecs(med_float  *val_med,
                             ecs_int_t   nbr_val)
{
  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(val_med != NULL);
  assert(nbr_val != 0);

  if (sizeof(med_float) == sizeof(ecs_coord_t)) {

    /* Les réels de type `med_float' et `ecs_coord_t' */
    /* sont codés sur le même nombre d'octets */

    /* Aucune conversion n'est nécessaire */

    return (ecs_coord_t *)val_med;
  }
  else {

    ecs_int_t   ival;
    ecs_coord_t    * val_ecs;

    /* Les réels de type `med_float' et `ecs_coord_t' */
    /* ne sont pas codés sur le même nombre d'octets  */

    /* On effectue la conversion de type pour chaque valeur */

    ECS_MALLOC(val_ecs, nbr_val, ecs_coord_t);

    for (ival = 0; ival < nbr_val; ival++)
      val_ecs[ival] = (ecs_coord_t)val_med[ival];

    ECS_FREE(val_med);

    return val_ecs;
  }
}

/*----------------------------------------------------------------------------
 *                        Lecture des noeuds
 *----------------------------------------------------------------------------*/

static void
ecs_loc_pre_med__lit_noeud(ecs_maillage_t   *maillage,
                           const ecs_med_t  *fic_maillage,
                           const char       *nom_maillage,
                           int               dim_e)
{
  /* Declarations des variables de stockage        */
  /* avant transfert dans la structure du maillage */
  /*-----------------------------------------------*/

  ecs_int_t     nbr_som;
  ecs_coord_t  *som_val_coord;

  /* Declarations des variables pour MED */
  /*-------------------------------------*/

  med_int      mdim_med;

  med_int      nbr_noe_med;
  med_int     *fam_noe_med;
  med_int     *num_noe_med;

  med_float   *coord_med;

  med_err      ret_med = 0;

  med_bool     changement = MED_FALSE;
  med_bool     transformation = MED_FALSE;

  char         nom_maillage_med[MED_NAME_SIZE + 1];

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  strcpy(nom_maillage_med, nom_maillage);

  /* On recupere le nombre de noeuds */

  nbr_noe_med = MEDmeshnEntity(fic_maillage->fid,
                               nom_maillage_med,
                               MED_NO_DT,
                               MED_NO_IT,
                               MED_NODE,
                               MED_NONE,
                               MED_COORDINATE,
                               MED_NODAL,
                               &changement,
                               &transformation);

  if (nbr_noe_med < 0)
    ret_med = -1;

  if (ret_med != 0)
    ecs_error(__FILE__, __LINE__, 0,
              _("MED: error reading file \"%s\".\n"
                "Number of vertices of entity read: \"%d\"."),
              fic_maillage->nom_fic, (int)nbr_noe_med);

  /*------------------------------------*/
  /* Lecture des coordonnées des noeuds */
  /*------------------------------------*/

  mdim_med = (med_int)dim_e;

  ECS_MALLOC(coord_med,   nbr_noe_med * mdim_med,            med_float);
  ECS_MALLOC(fam_noe_med, nbr_noe_med,                       med_int);
  ECS_MALLOC(num_noe_med, nbr_noe_med,                       med_int);

  ret_med = MEDmeshNodeCoordinateRd(fic_maillage->fid,
                                    nom_maillage_med,
                                    MED_NO_DT,
                                    MED_NO_IT,
                                    MED_FULL_INTERLACE,
                                    coord_med);

  if (ret_med != 0)
    ecs_error(__FILE__, __LINE__, 0,
              _("MED: error reading file \"%s\".\n"
                "Error reading coordinates."),
              fic_maillage->nom_fic);

  ECS_FREE(fam_noe_med);
  ECS_FREE(num_noe_med);

  /* Transformation du tableau des coord au type `med_float'   */
  /*             en un tableau des coord au type `ecs_coord_t' */
  /*-----------------------------------------------------------*/

  nbr_som = (ecs_int_t)nbr_noe_med;

  som_val_coord = ecs_loc_convert_real_med_ecs(coord_med,
                                               nbr_som * dim_e);

  if (dim_e != 3) {

    ecs_int_t isom;

    ECS_REALLOC(som_val_coord, nbr_som * 3, ecs_coord_t);

    if (dim_e == 1) {

      for (isom = nbr_som-1; isom >= 0; isom--) {
        som_val_coord[isom*3    ] = som_val_coord[isom];
        som_val_coord[isom*3 + 1] = 0.0;
        som_val_coord[isom*3 + 2] = 0.0;
      }

    }
    else if (dim_e == 2) {

      for (isom = nbr_som-1; isom >= 0; isom--) {
        som_val_coord[isom*3    ] = som_val_coord[isom*2    ];
        som_val_coord[isom*3 + 1] = som_val_coord[isom*2 + 1];
        som_val_coord[isom*3 + 2] = 0.0;
      }

    }
  }

  /* Transfert des valeurs lues dans la structure d'entite de maillage */
  /*===================================================================*/

  ecs_maillage_pre__cree_som(maillage,
                             nbr_som,
                             som_val_coord);
}

/*----------------------------------------------------------------------------
 *                        Lecture des éléments
 *----------------------------------------------------------------------------*/

static void
ecs_loc_pre_med__lit_maille(ecs_maillage_t   *maillage,
                            const ecs_med_t  *fic_maillage,
                            const char       *nom_maillage)
{
  ecs_entmail_t   entmail_e;

  ecs_int_t     cpt_som;
  ecs_int_t     ient;
  ecs_int_t     ielt;
  ecs_int_t     isom;
  ecs_int_t     ityp;
  ecs_int_t     nbr_elt;
  ecs_int_t     nbr_som_elt;
  ecs_int_t     nbr_som_fac;
  ecs_int_t     pos_elt;
  ecs_int_t     renum_som;
  ecs_int_t     taille;
  ecs_int_t     typ_geo_ecs;

  ecs_int_t    *elt_val_som;

  /* Déclarations des variables de stockage        */
  /* avant transfert dans la structure du maillage */
  /*-----------------------------------------------*/

  size_t      cpt_elt_ent[ECS_N_ENTMAIL];        /* Nbr. elts par entite */
  ecs_size_t *elt_pos_som_ent[ECS_N_ENTMAIL];    /* Positions sommets    */
  ecs_int_t  *elt_val_som_ent[ECS_N_ENTMAIL];    /* Numeros des sommets  */
  int        *elt_val_fam_ent[ECS_N_ENTMAIL];    /* Familles elements    */

  /* Déclarations des variables pour MED */
  /*-------------------------------------*/

  char               nom_maillage_med[MED_NAME_SIZE + 1] = "";
  char               nom_equiv[MED_NAME_SIZE+1] = "";

  med_geometry_type  typ_geo_med;
  med_entity_type    typ_ent_med;

  med_err            ret_med = 0;

  med_int            nbr_equiv = 0;
  med_int            ind_equiv = -1;
  med_int            edim_med;
  med_int            nbr_ele_med;
  med_int            taille_med;

  med_int          * connect_med = NULL;
  med_int          * fam_ele_med = NULL;

  med_bool           changement = MED_FALSE;
  med_bool           transformation = MED_FALSE;
  med_int            nstep = 0;
  med_int            nocstpcorrespondence = 0;
  med_data_type      data_type = MED_CONNECTIVITY;

  const int          n_typ_geo_ignore = 3;
  med_geometry_type  typ_geo_med_ignore[] = {MED_POINT1,
                                             MED_SEG2,
                                             MED_SEG3};
  const char        *nom_typ_ignore[] = {"point1", "seg2", "seg3"};

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /*====================================================*/
  /* Initialisations et allocations des tableaux locaux */
  /*====================================================*/

  /* Attention au decalage de `1' !!!         */
  /* On n'alloue pas les tableaux locaux pour */
  /* `ECS_ENTMAIL_DEB = ECS_ENTMAIL_SOM'      */

  for (ient = ECS_ENTMAIL_FAC; ient < ECS_N_ENTMAIL; ient++) {

    cpt_elt_ent        [ient] = 0;

    elt_pos_som_ent    [ient] = NULL;
    elt_val_som_ent    [ient] = NULL;
    elt_val_fam_ent    [ient] = NULL;

  }

  strcpy(nom_maillage_med, nom_maillage);

  /*------------------------------------------------*/
  /* Lecture de la connectivité nodale des elements */
  /*------------------------------------------------*/

  for (ityp = 0; ityp < n_typ_geo_ignore; ityp++) {

    /* On regarde si le maillage MED comporte des éléments ignorés */
    /*-------------------------------------------------------------*/

    nbr_ele_med = MEDmeshnEntity(fic_maillage->fid,
                                 nom_maillage_med,
                                 MED_NO_DT,
                                 MED_NO_IT,
                                 MED_NODE,
                                 typ_geo_med_ignore[ityp],
                                 MED_CONNECTIVITY,
                                 MED_NODAL,
                                 &changement,
                                 &transformation);

    if (nbr_ele_med < 0)
      ret_med = -1;

    else if (nbr_ele_med == 0) {

      nbr_ele_med = MEDmeshnEntity(fic_maillage->fid,
                                   nom_maillage_med,
                                   MED_NO_DT,
                                   MED_NO_IT,
                                   MED_CELL,
                                   typ_geo_med_ignore[ityp],
                                   MED_CONNECTIVITY,
                                   MED_NODAL,
                                   &changement,
                                   &transformation);

      if (nbr_ele_med < 0)
        ret_med = -1;

    }

    if (ret_med != 0)
      ecs_error(__FILE__, __LINE__, 0,
                _("MED: error reading file \"%s\".\n"
                  "Number of elements of entity read: \"%d\"."),
                fic_maillage->nom_fic, (int)nbr_ele_med);

    if (nbr_ele_med != 0) {

      ecs_warn();
      printf(_("The MED mesh contains %d elements of type %s\n"
               "which are ignored by the Preprocessor.\n\n"),
             (int)nbr_ele_med, nom_typ_ignore[ityp]);

    }
  }

  /* Vérification du nombre d'équivalences
     (pour connectivité faces non conformes éventuelle) */

  nbr_equiv = MEDnEquivalence(fic_maillage->fid,
                              nom_maillage_med);

  if (nbr_equiv < 0)
    ecs_error(__FILE__, __LINE__, 0,
              _("MED: error reading file \"%s\".\n"
                "Error reading equivalence information."),
              fic_maillage->nom_fic);

  else if (nbr_equiv > 0) {

    size_t  ind;
    char nom_equiv_cmp[MED_NAME_SIZE+1];
    char desc_equiv[MED_COMMENT_SIZE+1];

    for (ind_equiv = 0; ind_equiv < nbr_equiv; ind_equiv++) {

      ret_med = MEDequivalenceInfo(fic_maillage->fid,
                                   nom_maillage_med,
                                   ind_equiv + 1,
                                   nom_equiv,
                                   desc_equiv,
                                   &nstep,
                                   &nocstpcorrespondence);

      if (ret_med < 0)
        ecs_error(__FILE__, __LINE__, 0,
                  _("MED: error reading file \"%s\".\n"
                    "Error reading equivalence information."),
                  fic_maillage->nom_fic);

      nom_equiv[MED_NAME_SIZE] = '\0';

      for (ind = 0; ind < MED_NAME_SIZE && nom_equiv[ind] != '\0'; ind++)
        nom_equiv_cmp[ind] = tolower(nom_equiv[ind]);
      nom_equiv_cmp[strlen(nom_equiv)] = '\0';

      if (   strcmp(nom_equiv_cmp, "face_connectivity") == 0
          || strcmp(nom_equiv_cmp, "face connectivity") == 0)
        break;

    }

    if (ind_equiv >= nbr_equiv) {
      nom_equiv[0] = '\0';
      ind_equiv = -1;
    }

  }

  /* Boucle sur tous les types géométriques MED référencés dans l'Enveloppe */
  /*------------------------------------------------------------------------*/

  for (ityp = 0; ityp < ECS_MED_NBR_TYP_ELT; ityp++) {

    taille_med = 0;
    typ_geo_med = ecs_fic_med_init_elt_liste_c[ityp].med_typ;

    /*
     * On essaiera de lire en priorité les mailles, les arêtes ou faces
     * correspondant à une connectivité descendante étant lues en second choix.
     */

    edim_med = typ_geo_med / 100;

    typ_ent_med = MED_CELL;

    nbr_ele_med = MEDmeshnEntity(fic_maillage->fid,
                                 nom_maillage_med,
                                 MED_NO_DT,
                                 MED_NO_IT,
                                 typ_ent_med,
                                 typ_geo_med,
                                 MED_CONNECTIVITY,
                                 MED_NODAL,
                                 &changement,
                                 &transformation);

    if (nbr_ele_med > 0) { /* Special case for polygons and polyhedra */

      if (typ_geo_med == MED_POLYGON)
        data_type = MED_INDEX_NODE;
      else if (typ_geo_med == MED_POLYHEDRON)
        data_type = MED_INDEX_FACE;

      if (data_type != MED_CONNECTIVITY)
        nbr_ele_med = MEDmeshnEntity(fic_maillage->fid,
                                     nom_maillage_med,
                                     MED_NO_DT,
                                     MED_NO_IT,
                                     typ_ent_med,
                                     typ_geo_med,
                                     data_type,
                                     MED_NODAL,
                                     &changement,
                                     &transformation) - 1;

    }

    if (nbr_ele_med < 0)
      ret_med = -1;

    else if (nbr_ele_med == 0 && edim_med <= 2) {

      if (edim_med == 1)
        typ_ent_med = MED_DESCENDING_EDGE;
      else if (edim_med == 2)
        typ_ent_med = MED_DESCENDING_FACE;

      nbr_ele_med = MEDmeshnEntity(fic_maillage->fid,
                                   nom_maillage_med,
                                   MED_NO_DT,
                                   MED_NO_IT,
                                   typ_ent_med,
                                   typ_geo_med,
                                   MED_CONNECTIVITY,
                                   MED_NODAL,
                                   &changement,
                                   &transformation);

      if (nbr_ele_med < 0)
        ret_med = -1;

    }

    if (ret_med != 0)
      ecs_error(__FILE__, __LINE__, 0,
                _("MED: error reading file \"%s\".\n"
                  "Number of elements of entity read: \"%d\"."),
                fic_maillage->nom_fic, (int)nbr_ele_med);

    if (nbr_ele_med != 0) {

      nbr_elt = (ecs_int_t)nbr_ele_med;

      /* Type géometrique */
      /*------------------*/

      typ_geo_ecs = ecs_fic_med_init_elt_liste_c[ityp].ecs_typ;

      /* Identification de l'entité concernée */
      /*--------------------------------------*/

      entmail_e = ecs_maillage_pre__ret_typ_geo(typ_geo_ecs);

      ECS_REALLOC(elt_pos_som_ent[entmail_e],
                  cpt_elt_ent[entmail_e] + 1 + nbr_elt, ecs_size_t);

      if (cpt_elt_ent[entmail_e] == 0) {
        /* On est au 1er tour */
        elt_pos_som_ent[entmail_e][0] = 1;
      }

      ECS_REALLOC(elt_val_fam_ent[entmail_e],
                  cpt_elt_ent[entmail_e] + nbr_elt, int);

      ECS_MALLOC(fam_ele_med, nbr_ele_med, med_int);

      /* Traitement des éléments "classiques" */
      /*--------------------------------------*/

      if (typ_geo_med != MED_POLYGON && typ_geo_med != MED_POLYHEDRON) {

        nbr_som_elt = ecs_fic_elt_typ_liste_c[typ_geo_ecs].nbr_som;

        taille = elt_pos_som_ent[entmail_e][cpt_elt_ent[entmail_e]] - 1;

        ECS_REALLOC(elt_val_som_ent[entmail_e],
                    taille + nbr_elt * nbr_som_elt, ecs_int_t);


        /* Convention sur la taille des mailles */

        edim_med = typ_geo_med / 100;
        taille_med = typ_geo_med % 100;

        taille = (ecs_int_t)nbr_ele_med * (ecs_int_t)taille_med;
        ECS_MALLOC(connect_med, taille, med_int);

        ret_med = MEDmeshElementConnectivityRd(fic_maillage->fid,
                                               nom_maillage_med,
                                               MED_NO_DT,
                                               MED_NO_IT,
                                               typ_ent_med,
                                               typ_geo_med,
                                               MED_NODAL,
                                               MED_FULL_INTERLACE,
                                               connect_med);

        if (ret_med != 0)
          ecs_error(__FILE__, __LINE__, 0,
                    _("MED: error reading file \"%s\".\n"
                      "Error reading connectivity."),
                    fic_maillage->nom_fic);

        /* Les références aux noeuds milieux des éléments     */
        /*  quadratiques ne sont pas conservées               */

        elt_val_som = ecs_loc_convert_connect_med_ecs(connect_med,
                                                      taille,
                                                      taille_med,
                                                      nbr_som_elt);

        /* Remplissage de la connectivité */

        cpt_som = 0;

        for (ielt = 0; ielt < nbr_elt; ielt++) {

          pos_elt = elt_pos_som_ent[entmail_e][cpt_elt_ent[entmail_e]
                                               + ielt];

          elt_pos_som_ent[entmail_e][cpt_elt_ent[entmail_e] + ielt + 1]
            = pos_elt + nbr_som_elt;

          for (isom = 0; isom < nbr_som_elt; isom++) {

            renum_som = ecs_fic_med_init_elt_liste_c[ityp].num_som[isom] - 1;
            elt_val_som_ent[entmail_e][pos_elt - 1 + isom]
              = elt_val_som[cpt_som + renum_som];

          }

          cpt_som += nbr_som_elt;
        }

        /* Libération connectivité temporaire */

        ECS_FREE(elt_val_som);
      }

      /* Traitement des polygones */
      /*--------------------------*/

      else if (typ_geo_med == MED_POLYGON) {

        ecs_int_t    ival;
        ecs_int_t    nbr_val_elt;
        med_int    * index_med = NULL;

        /* Taille du tableau des connectivites */

        taille_med = MEDmeshnEntity(fic_maillage->fid,
                                    nom_maillage_med,
                                    MED_NO_DT,
                                    MED_NO_IT,
                                    typ_ent_med,
                                    MED_POLYGON,
                                    MED_CONNECTIVITY,
                                    MED_NODAL,
                                    &changement,
                                    &transformation);

        if (taille_med < 0)
          ret_med = -1;

        if (ret_med != 0)
          ecs_error(__FILE__, __LINE__, 0,
                    _("MED: error reading file \"%s\".\n"
                      "(polygons information)."),
                    fic_maillage->nom_fic);

        taille = elt_pos_som_ent[entmail_e][cpt_elt_ent[entmail_e]] - 1;

        ECS_REALLOC(elt_val_som_ent[entmail_e],
                    taille + (ecs_int_t)taille_med, ecs_int_t);

        ECS_MALLOC(index_med, (ecs_int_t)(nbr_ele_med + 1), med_int);
        ECS_MALLOC(connect_med, (ecs_int_t)taille_med, med_int);

        /* Lecture de la connectivité des polygones */

        ret_med = MEDmeshPolygonRd(fic_maillage->fid,
                                   nom_maillage_med,
                                   MED_NO_DT,
                                   MED_NO_IT,
                                   typ_ent_med,
                                   MED_NODAL,
                                   index_med,
                                   connect_med);

        if (ret_med != 0)
          ecs_error(__FILE__, __LINE__, 0,
                    _("MED: error reading file \"%s\".\n"
                      "(polygons connectivity)."),
                    fic_maillage->nom_fic);

        /* Remplissage de la connectivité */

        pos_elt = elt_pos_som_ent[entmail_e][cpt_elt_ent[entmail_e]];

        for (ielt = 0; ielt < nbr_elt; ielt++)
          elt_pos_som_ent[entmail_e][cpt_elt_ent[entmail_e] + ielt + 1]
            = pos_elt + index_med[ielt + 1] - index_med[0];

        nbr_val_elt = index_med[nbr_elt] - index_med[0];

        for (ival = 0; ival < nbr_val_elt; ival++)
          elt_val_som_ent[entmail_e][pos_elt - 1 + ival]
            = connect_med[ival];

        /* Libération connectivité temporaire */

        ECS_FREE(index_med);
        ECS_FREE(connect_med);

      }

      /* Traitement des polyèdres */
      /*--------------------------*/

      else if (typ_geo_med == MED_POLYHEDRON) {

        ecs_int_t    ifac, ival;
        ecs_int_t    num_som_prec, num_som_deb;
        ecs_int_t    cpt_som_poly_dup = 0;
        ecs_int_t    cpt_fac_poly_dgn = 0;
        med_int      n_index_f_med = 0;
        med_int    * index_med = NULL;
        med_int    * index_f_med = NULL;

        /* Taille du tableau des connectivites */

        ret_med = 0;

        n_index_f_med = MEDmeshnEntity(fic_maillage->fid,
                                       nom_maillage_med,
                                       MED_NO_DT,
                                       MED_NO_IT,
                                       typ_ent_med,
                                       MED_POLYHEDRON,
                                       MED_INDEX_NODE,
                                       MED_NODAL,
                                       &changement,
                                       &transformation);

        if (n_index_f_med < 0)
          ret_med = -1;

        if (ret_med >= 0) {

          taille_med = MEDmeshnEntity(fic_maillage->fid,
                                      nom_maillage_med,
                                      MED_NO_DT,
                                      MED_NO_IT,
                                      typ_ent_med,
                                      MED_POLYHEDRON,
                                      MED_CONNECTIVITY,
                                      MED_NODAL,
                                      &changement,
                                      &transformation);

          if (taille_med < 0)
            ret_med = -1;

        }

        if (ret_med != 0)
          ecs_error(__FILE__, __LINE__, 0,
                    _("MED: error reading file \"%s\".\n"
                      "(polyhedra information)."),
                    fic_maillage->nom_fic);

        /* On ajoute le nombre de faces (n_index_f_med - 1) à la taille du
           tableau destiné à recevoir la connectivité afin de faire apparaître
           une seconde fois le numéro du premier sommet à la fin de la liste
           des sommets de chaque face dans la définition des cellules (pour
           repérer la fin de définition de chaque face d'une cellule). */

        taille = elt_pos_som_ent[entmail_e][cpt_elt_ent[entmail_e]] - 1;

        ECS_REALLOC(elt_val_som_ent[entmail_e],
                    taille + (ecs_int_t)(taille_med + n_index_f_med - 1),
                    ecs_int_t);

        ECS_MALLOC(index_med, (ecs_int_t)(nbr_ele_med + 1), med_int);
        ECS_MALLOC(index_f_med, (ecs_int_t)(n_index_f_med), med_int);
        ECS_MALLOC(connect_med, (ecs_int_t)taille_med, med_int);

        /* Lecture de la connectivité des polyèdres */

        ret_med = MEDmeshPolyhedronRd(fic_maillage->fid,
                                      nom_maillage_med,
                                      MED_NO_DT,
                                      MED_NO_IT,
                                      typ_ent_med,
                                      MED_NODAL,
                                      index_med,
                                      index_f_med,
                                      connect_med);

        if (ret_med != 0)
          ecs_error(__FILE__, __LINE__, 0,
                    _("MED: error reading file \"%s\".\n"
                      "(polyhedra connectivity)."),
                    fic_maillage->nom_fic);

        if (index_med[nbr_ele_med] - index_med[0] + 1 != n_index_f_med)
          ecs_error
            (__FILE__, __LINE__, 0,
             _("MED: inconsistency in polyhedra connectivity;\n"
               "number of polyhedra: %d\n"
               "the cells->faces end index (%d) should be equal to\n"
               "the size of the faces->vertices index array (%d)\n"),
             (int) nbr_ele_med,
             (int)(index_med[nbr_ele_med] - index_med[0] + 1),
             (int) n_index_f_med);

        /* Remplissage de la connectivité */

        pos_elt = elt_pos_som_ent[entmail_e][cpt_elt_ent[entmail_e]];

        ival = pos_elt - 1;

        for (ielt = 0; ielt < nbr_elt; ielt++) {

          for (ifac = index_med[ielt    ] - index_med[0];
               ifac < index_med[ielt + 1] - index_med[0];
               ifac++) {

            /* Premier sommet de la face */

            isom = index_f_med[ifac    ] - index_f_med[0];

            elt_val_som_ent[entmail_e][ival++] = connect_med[isom];

            num_som_prec = connect_med[isom];
            num_som_deb  = connect_med[isom];

            nbr_som_fac = 1;

            /* Autres sommets de la face */

            for (isom = index_f_med[ifac    ] - index_f_med[0] + 1;
                 isom < index_f_med[ifac + 1] - index_f_med[0];
                 isom++) {

              if (   connect_med[isom] != num_som_prec
                  && connect_med[isom] != num_som_deb) {
                elt_val_som_ent[entmail_e][ival++] = connect_med[isom];
                num_som_prec = connect_med[isom];
                nbr_som_fac += 1;
              }
              else
                cpt_som_poly_dup += 1;

            }

            /* Si la face est dégénérée (1 ou 2 sommets), on la supprime */

            if (nbr_som_fac < 3) {

              ival -= nbr_som_fac;
              cpt_fac_poly_dgn += 1;

            }
            else {

              /* On rajoute une seconde référence au premier sommet de
                 chaque face en fin de liste pour "fermer" cette face
                 (convention polyèdre Enveloppe Code_Saturne pour repérer
                 les faces d'un polyèdre en connectivité nodale) */

              isom = index_f_med[ifac    ] - index_f_med[0];

              elt_val_som_ent[entmail_e][ival++] = connect_med[isom];

            }

          }

          elt_pos_som_ent[entmail_e][cpt_elt_ent[entmail_e] + ielt + 1]
            = ival + 1;

        }

        assert((size_t)(taille_med + n_index_f_med - 1)
               >= (size_t)((  elt_pos_som_ent[entmail_e]
                                             [cpt_elt_ent[entmail_e] + ielt]
                            - elt_pos_som_ent[entmail_e]
                                             [cpt_elt_ent[entmail_e]])
                           + cpt_som_poly_dup + cpt_fac_poly_dgn));

        if (cpt_som_poly_dup > 0) {

          ecs_warn();
          printf(_("While reading file \"%s\",\n"
                   "%d repeated references to the same vertices\n"
                   "and %d degenerate faces were encountered\n"
                   "in the definition of %d polyhedra.\n"),
                 fic_maillage->nom_fic, (int) cpt_som_poly_dup,
                 (int) cpt_fac_poly_dgn, (int) nbr_elt);

          ECS_REALLOC (elt_val_som_ent[entmail_e], ival + 1, ecs_int_t);

        }

        /* Libération connectivité temporaire */

        ECS_FREE(index_med);
        ECS_FREE(index_f_med);
        ECS_FREE(connect_med);

      }

      /* Lecture des numéros de familles */

      ret_med = MEDmeshEntityFamilyNumberRd(fic_maillage->fid,
                                            nom_maillage_med,
                                            MED_NO_DT,
                                            MED_NO_IT,
                                            typ_ent_med,
                                            typ_geo_med,
                                            fam_ele_med);

      /* Convention MED : si pas de familles, numéros = 0 */

      if (ret_med < 0) {
        for (ielt = 0; ielt < nbr_elt; ielt++)
          fam_ele_med[ielt] = 0;
        ret_med = 0;
      }
      else if (ret_med > 0) {
        ecs_warn();
        printf(_("MED: erreur reading file \"%s\".\n"
                 "(families)."),
               fic_maillage->nom_fic);
      }

      /* Changement de signe pour les familles des éléments */

      for (ielt = 0; ielt < nbr_elt; ielt++) {

        elt_val_fam_ent[entmail_e][cpt_elt_ent[entmail_e] + ielt]
          = -fam_ele_med[ielt];

      }

      ECS_FREE(fam_ele_med);

      /* Lecture éventuelle des connectivités */

      if (edim_med == 2 && ind_equiv >= 0) {

        size_t  n_corres_old = maillage->n_connect_couples[ECS_ENTMAIL_FAC];
        med_int nbr_corres;
        med_int *corres = NULL;
        ecs_int_t *copy_fac_connect = NULL;

        ret_med = MEDequivalenceCorrespondenceSize(fic_maillage->fid,
                                                   nom_maillage_med,
                                                   nom_equiv,
                                                   MED_NO_DT,
                                                   MED_NO_IT,
                                                   typ_ent_med,
                                                   typ_geo_med,
                                                   &nbr_corres);

        if (ret_med < 0)
          ecs_error(__FILE__, __LINE__, 0,
                    _("MED: error reading file \"%s\".\n"
                      "Error reading equivalence information."),
                    fic_maillage->nom_fic);

        else if (nbr_corres > 0) {

          maillage->n_connect_couples[ECS_ENTMAIL_FAC] += nbr_corres;

          ECS_REALLOC(maillage->connect_couples[ECS_ENTMAIL_FAC],
                      maillage->n_connect_couples[ECS_ENTMAIL_FAC] * 2,
                      ecs_int_t);

          ECS_MALLOC(corres, nbr_corres*2, med_int);

          ret_med = MEDequivalenceCorrespondenceRd(fic_maillage->fid,
                                                   nom_maillage_med,
                                                   nom_equiv,
                                                   MED_NO_DT,
                                                   MED_NO_IT,
                                                   typ_ent_med,
                                                   typ_geo_med,
                                                   corres);

        }

        /* Add new connected faces, shifting numbering in case of
           previously read element types. */

        copy_fac_connect
          = maillage->connect_couples[ECS_ENTMAIL_FAC] + (n_corres_old*2);

        for (ielt = 0; ielt < nbr_corres*2; ielt++)
          copy_fac_connect[ielt] = corres[ielt] + cpt_elt_ent[entmail_e];

        ECS_FREE(corres);
      }

      cpt_elt_ent[entmail_e] += nbr_elt;

    } /* Fin : si le nombre d'éléments de ce type MED n'est pas nul */

  } /* Fin : boucle sur les types d'elements MED supportés */

  /* Transfert des valeurs lues dans les structures d'entité de maillage */
  /*=====================================================================*/

  ecs_maillage_pre__cree_elt(maillage,
                             cpt_elt_ent,
                             elt_pos_som_ent,
                             elt_val_som_ent,
                             elt_val_fam_ent,
                             NULL,
                             NULL,
                             NULL,
                             NULL);
}

/*----------------------------------------------------------------------------
 *                        Lecture des familles
 *----------------------------------------------------------------------------*/

static ecs_famille_t **
ecs_loc_pre_med__lit_famille(const ecs_med_t  *fic_maillage,
                             const char       *nom_maillage)
{
  char                *nom;

  ecs_int_t            deb_pos_sans_blanc;
  ecs_int_t            fin_pos_sans_blanc;
  ecs_int_t            icar;
  ecs_int_t            ide;
  ecs_int_t            ifam;
  ecs_int_t            ifam_ent;
  ecs_int_t            nbr_car;
  ecs_int_t            nbr_famille_ent[ECS_N_ENTMAIL];
  ecs_int_t            num_fam;
  ecs_int_t            num_fam_ent;
  ecs_descr_t         *descr;
  ecs_descr_t         *descr_tete;

  ecs_famille_t       *famille;
  ecs_famille_t     **vect_famille_tete;
  ecs_famille_t     **liste_famille_ent[ECS_N_ENTMAIL];

  /* Declarations des variables pour MED */
  /*-------------------------------------*/

  char       nom_maillage_med[MED_NAME_SIZE + 1];
  char      *att_des_med;
  char      *grp_des_med;
  char       nom_fam_med[MED_NAME_SIZE + 1];
  char       un_att_des_med[MED_COMMENT_SIZE + 1];
  char       un_grp_des_med[MED_LNAME_SIZE + 1];

  med_err    ret_med = 0;

  med_int    iatt_med;
  med_int    ifam_med;
  med_int    igrp_med;
  med_int    nbr_fam_med;
  med_int    nbr_att_med;
  med_int    nbr_grp_med;
  med_int    num_fam_med;
  med_int   *att_ide_med;
  med_int   *att_val_med;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  strcpy(nom_maillage_med, nom_maillage);

  /* On récupere le nombre de familles */

  nbr_fam_med = MEDnFamily(fic_maillage->fid,
                           nom_maillage_med);

  if (nbr_fam_med < 0)
    ret_med = -1;

  if (ret_med != 0)
    ecs_error(__FILE__, __LINE__, 0,
              _("MED: error reading file \"%s\".\n"
                "Number of families read: \"%d\""),
              fic_maillage->nom_fic, (int)nbr_fam_med);

  for (ifam_ent = 0; ifam_ent < ECS_N_ENTMAIL; ifam_ent++) {

    nbr_famille_ent[ifam_ent] = 0;
    ECS_MALLOC(liste_famille_ent[ifam_ent],
               (ecs_int_t)nbr_fam_med, ecs_famille_t *);

  }

  ECS_MALLOC(vect_famille_tete, ECS_N_ENTMAIL, ecs_famille_t *);
  for (ifam_ent = 0; ifam_ent < ECS_N_ENTMAIL; ifam_ent++)
    vect_famille_tete[ifam_ent] = NULL;

  /*----------------------*/
  /* Lecture des familles */
  /*----------------------*/

  for (ifam_med = 0; ifam_med < nbr_fam_med; ifam_med++) {

    /*------------------------------------*/
    /* Récupération du nombre d'attributs */
    /*------------------------------------*/

    if (fic_maillage->version[0] == 2)
      nbr_att_med = MEDnFamily23Attribute(fic_maillage->fid,
                                          nom_maillage_med,
                                          ifam_med + 1);
    else
      nbr_att_med = 0;

    if (nbr_att_med < 0)
      ret_med = -1;

    if (ret_med != 0)
      ecs_error(__FILE__, __LINE__, 0,
                _("MED: error reading file \"%s\".\n"
                  "Number of attributes read: \"%d\""),
                fic_maillage->nom_fic, (int)nbr_att_med);

    /*-----------------------------------*/
    /* Récupération du nombre de groupes */
    /*-----------------------------------*/

    nbr_grp_med = MEDnFamilyGroup(fic_maillage->fid,
                                  nom_maillage_med,
                                  ifam_med + 1);

    if (nbr_grp_med < 0)
      ret_med = -1;

    if (ret_med != 0)
      ecs_error(__FILE__, __LINE__, 0,
                _("MED: error reading file \"%s\".\n"
                  "Number of groups read: \"%d\""),
                fic_maillage->nom_fic, (int)nbr_grp_med);

    /*--------------------------------------*/
    /* Lecture des attributs et des groupes */
    /*--------------------------------------*/

    if (nbr_att_med != 0) {

      ECS_MALLOC(att_ide_med, nbr_att_med, med_int);
      ECS_MALLOC(att_val_med, nbr_att_med, med_int);
      ECS_MALLOC(att_des_med, nbr_att_med * MED_COMMENT_SIZE + 1, char);

    }
    else {

      att_ide_med = NULL;
      att_val_med = NULL;
      att_des_med = NULL;

    }

    ECS_MALLOC(grp_des_med, nbr_grp_med * MED_LNAME_SIZE + 1, char);

    if (fic_maillage->version[0] == 2)
      ret_med = MEDfamily23Info(fic_maillage->fid,
                                nom_maillage_med,
                                ifam_med + 1,
                                nom_fam_med,
                                att_ide_med,
                                att_val_med,
                                att_des_med,
                                &num_fam_med,
                                grp_des_med);

    else
      ret_med = MEDfamilyInfo(fic_maillage->fid,
                              nom_maillage_med,
                              ifam_med + 1,
                              nom_fam_med,
                              &num_fam_med,
                              grp_des_med);

    if (num_fam_med != 0) {

      descr_tete = NULL;

      if (nbr_att_med == 0 && nbr_grp_med == 0)
        printf(_("  Family %2d is described by no attribute or group.\n"),
               (int)(ECS_ABS(num_fam_med)));

      for (iatt_med = 0; iatt_med < nbr_att_med; iatt_med++) {

        /* Récupération de la valeur entière du descripteur */
        /*--------------------------------------------------*/

        ide =  (ecs_int_t)att_val_med[iatt_med];

        /* Recupération du descripteur */
        /*-----------------------------*/

        strncpy(un_att_des_med,
                att_des_med + iatt_med * MED_COMMENT_SIZE, MED_COMMENT_SIZE);
        un_att_des_med[MED_COMMENT_SIZE] = '\0';

        /* On regarde si la chaîne ne contient pas que des blancs */
        /* On ne garde que la partie de la chaîne                 */
        /*  qui n'a pas de blancs aux extrémités                  */

        deb_pos_sans_blanc = 0;
        while (*(un_att_des_med + deb_pos_sans_blanc) != '\0' &&
               *(un_att_des_med + deb_pos_sans_blanc) == ' ')
          deb_pos_sans_blanc++;

        for (fin_pos_sans_blanc = deb_pos_sans_blanc;
             (   fin_pos_sans_blanc < MED_COMMENT_SIZE - 1
               && *(un_att_des_med + fin_pos_sans_blanc) != '\0');
             fin_pos_sans_blanc++);

        if (   fin_pos_sans_blanc > deb_pos_sans_blanc
            && deb_pos_sans_blanc < MED_COMMENT_SIZE) {

          /* La chaîne ne contient pas que des blancs */

          while (fin_pos_sans_blanc                     != 0    &&
                 *(un_att_des_med + fin_pos_sans_blanc) == ' ')
            fin_pos_sans_blanc--;

          nbr_car = fin_pos_sans_blanc - deb_pos_sans_blanc + 1;
          ECS_MALLOC(nom, nbr_car + 1, char);
          for (icar = 0; icar < nbr_car; icar++)
            *(nom + icar) = *(un_att_des_med + deb_pos_sans_blanc + icar);
          *(nom + nbr_car) = '\0';

        }
        else {
          nom   = NULL;
        }

        descr = ecs_descr__cree(ide, nom);

        ecs_descr_chaine__ajoute(&descr_tete, descr);

        if (nom != NULL)
          ECS_FREE(nom);

      } /* Fin : boucle sur les attributs de la famille */

      for (igrp_med = 0; igrp_med < nbr_grp_med; igrp_med++) {

        /* Récuperation de la chaîne de caracteres du descripteur */

        strncpy(un_grp_des_med,
                grp_des_med + igrp_med * MED_LNAME_SIZE, MED_LNAME_SIZE);
        un_grp_des_med[MED_LNAME_SIZE] = '\0';

        /* On regarde si la chaîne ne contient pas que des blancs */
        /* On ne garde que la partie de la chaîne                 */
        /*  qui n'a pas de blancs aux extrémités                  */

        deb_pos_sans_blanc = 0;
        while (*(un_grp_des_med + deb_pos_sans_blanc) != '\0' &&
               *(un_grp_des_med + deb_pos_sans_blanc) == ' ')
          deb_pos_sans_blanc++;

        /* On verifie que la chaîne ne contient pas que des blancs */
        assert(deb_pos_sans_blanc < MED_LNAME_SIZE + 1);

        fin_pos_sans_blanc = MED_LNAME_SIZE - 1;
        while (fin_pos_sans_blanc              != 0    &&
               *(un_grp_des_med + fin_pos_sans_blanc) == ' ')
          fin_pos_sans_blanc--;

        nbr_car = fin_pos_sans_blanc - deb_pos_sans_blanc + 1;
        ECS_MALLOC(nom, nbr_car + 1, char);
        for (icar = 0; icar < nbr_car; icar++)
          *(nom + icar) = *(un_grp_des_med + deb_pos_sans_blanc + icar);
        *(nom + nbr_car) = '\0';


        /* Pas de valeur entière associée */

        descr = ecs_descr__cree(ECS_DESCR_IDE_NUL, nom);

        ecs_descr_chaine__ajoute(&descr_tete, descr);

        if (nom != NULL)
          ECS_FREE(nom);

      } /* Fin : boucle sur les groupes de la famille */


      /* Détermination de l'entité concernée par la famille */

      if (num_fam_med < 0) {

        /* On accroche toutes les familles des éléments */
        /*  sur les familles des cellules               */

        num_fam_ent = ECS_ENTMAIL_CEL;

        num_fam = (ecs_int_t)num_fam_med;

        famille = ecs_famille__cree(num_fam, descr_tete);

        liste_famille_ent[num_fam_ent][nbr_famille_ent[num_fam_ent]] = famille;
        nbr_famille_ent[num_fam_ent]++;

      }
      else {

        num_fam_ent = ECS_ENTMAIL_NONE;

      }

    } /* Fin : si la famille n'est pas la famille par défaut de MED */

    if (nbr_att_med != 0) {

      ECS_FREE(att_ide_med);
      ECS_FREE(att_val_med);
      ECS_FREE(att_des_med);

    }
    ECS_FREE(grp_des_med);

  } /* Fin : boucle sur les familles lues */

  for (ifam_ent = 0; ifam_ent < ECS_N_ENTMAIL; ifam_ent++) {

    for (ifam = 0; ifam < nbr_famille_ent[ifam_ent]; ifam++) {

      ecs_famille_chaine__ajoute(&vect_famille_tete[ifam_ent],
                                 liste_famille_ent[ifam_ent][ifam]);
    }

  }

  for (ifam_ent = 0; ifam_ent < ECS_N_ENTMAIL; ifam_ent++)
    ECS_FREE(liste_famille_ent[ifam_ent]);

  return vect_famille_tete;
}

/*============================================================================
 *                             Fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Lecture d'un fichier au format MED
 *   et affectation des donnees dans la structure de maillage
 *----------------------------------------------------------------------------*/

ecs_maillage_t *
ecs_pre_med__lit_maillage(const char  *nom_fic_maillage,
                          int          num_maillage)
{
  char            *nom_maillage;
  int              ind;
  med_int          dim_e;
  ecs_med_t       *fic_maillage;
  ecs_famille_t  **vect_famille;

  /* Declarations des variables pour MED */
  /*-------------------------------------*/

  char           nom_maillage_med[MED_NAME_SIZE + 1];
  char           desc_maillage_med[MED_COMMENT_SIZE + 1];

  med_err        ret_med = 0;

  med_int        mdim_med;
  med_int        nbr_maillages_med;
  int            num_maillage_med = 1;
  med_mesh_type  type_maillage_med;

  char              dtunit[MED_LNAME_SIZE + 1];
  char              axisname[MED_SNAME_SIZE*3 + 1];
  char              axisunit[MED_SNAME_SIZE*3 + 1];
  med_sorting_type  sortingtype;
  med_int           nstep;
  med_axis_type     axistype;

  /* Création d'un maillage initialement vide (valeur de retour) */

  ecs_maillage_t  *maillage = ecs_maillage__cree_nodal();

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  printf(_("\n\n"
           "Reading mesh from file in MED (EDF/CEA) format\n"
           "----------------------\n"));

  printf(_("  Mesh file: %s\n\n\n"),
         nom_fic_maillage);

  /* Ouverture du fichier MED en lecture */
  /*-------------------------------------*/

  fic_maillage = ecs_pre_med__cree(nom_fic_maillage);

  /* Vérification et traitements selon le nombre de maillages */

  nbr_maillages_med = MEDnMesh(fic_maillage->fid);

  if (nbr_maillages_med < 0)
    ecs_error(__FILE__, __LINE__, 0,
              _("MED: error reading file \"%s\"."),
              nom_fic_maillage);

  if (nbr_maillages_med > 1) {
    printf(_("\n  The file contains multiple meshes:\n"));
    for (ind = 0; ind < nbr_maillages_med; ind++) {
      ret_med = MEDmeshInfo(fic_maillage->fid,
                            ind + 1,
                            nom_maillage_med,
                            &dim_e,
                            &mdim_med,
                            &type_maillage_med,
                            desc_maillage_med,
                            dtunit,
                            &sortingtype,
                            &nstep,
                            &axistype,
                            axisname,
                            axisunit);
      if (ret_med != 0)
        ecs_error(__FILE__, __LINE__, 0,
                  _("MED: error reading file \"%s\"."),
                  nom_fic_maillage);
      printf(_("    Mesh %2d: %s\n"), ind + 1, nom_maillage_med);
    }
    if (num_maillage == 0)
      printf(_("\n  No mesh was specified; the first is read.\n\n"));
    else if (num_maillage > 0)
      printf(_("\n  Mesh number %d was specified.\n\n"), num_maillage);
  }
  else if (nbr_maillages_med == 0)
    ecs_error(__FILE__, __LINE__, 0,
              _("No mesh in file.\n"));

  assert (num_maillage >= 0);
  if (num_maillage > nbr_maillages_med)
    ecs_error(__FILE__, __LINE__, 0,
              _("The specified mesh number (%d) is greater than\n"
                "the number of meshes defined (%d) in file\n%s.\n"),
              num_maillage, nbr_maillages_med, nom_fic_maillage);
  else
    num_maillage_med = ECS_MAX(1, num_maillage);

  ret_med = MEDmeshInfo(fic_maillage->fid,
                        num_maillage_med,
                        nom_maillage_med,
                        &dim_e,
                        &mdim_med,
                        &type_maillage_med,
                        desc_maillage_med,
                        dtunit,
                        &sortingtype,
                        &nstep,
                        &axistype,
                        axisname,
                        axisunit);

  if (ret_med != 0)
    ecs_error(__FILE__, __LINE__, 0,
              _("MED: error reading file \"%s\".\n"
                "Name of mesh read : \"%s\"\n"
                "Dimension read    : \"%d\""),
              nom_fic_maillage, nom_maillage_med, mdim_med);

  nom_maillage_med[MED_NAME_SIZE] = '\0';

  printf(_("  Mesh name: %s\n\n"), nom_maillage_med);

  assert((int)mdim_med == 2 || (int)mdim_med == 3);

  ECS_MALLOC(nom_maillage, strlen(nom_maillage_med) + 1, char);
  strcpy(nom_maillage, nom_maillage_med);

  /* Lecture des noeuds */
  /*--------------------*/

  ecs_loc_pre_med__lit_noeud(maillage,
                             fic_maillage,
                             nom_maillage,
                             dim_e);

  /* Lecture des elements */
  /*----------------------*/

  ecs_loc_pre_med__lit_maille(maillage,
                              fic_maillage,
                              nom_maillage);

  /* Lecture des familles */
  /*----------------------*/

  vect_famille = ecs_loc_pre_med__lit_famille(fic_maillage,
                                              nom_maillage);

  ecs_maillage__definit_famille(maillage,
                                vect_famille);

  ECS_FREE(vect_famille);
  ECS_FREE(nom_maillage);

  /* Fermeture du fichier de lecture du maillage */
  /*---------------------------------------------*/

  ecs_pre_med__detruit(fic_maillage);

  /* Transformation des familles en attributs "groupe" */
  /*---------------------------------------------------*/

  ecs_maillage__cree_attributs(maillage);

  return maillage;
}

/*----------------------------------------------------------------------------*/

#endif /* HAVE_MED */
