/*============================================================================
 *  Définition de la fonction
 *   de lecture d'un fichier de maillage au format CGNS
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
 * Visibilité
 *============================================================================*/


#include "cs_config.h"

#if defined(HAVE_CGNS)


/*----------------------------------------------------------------------------
 * Fichiers `include' librairie standard C
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>


/*----------------------------------------------------------------------------
 * Fichiers `include' publics  du  paquetage global "CGNS"
 *----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif

#include <cgnslib.h>

#ifdef __cplusplus
}
#endif


/*----------------------------------------------------------------------------
 * Fichiers `include' visibles du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/

#include "ecs_def.h"
#include "ecs_elt_typ_liste.h"
#include "ecs_mem.h"


/*----------------------------------------------------------------------------
 * Fichiers `include' visibles des paquetages visibles
 *----------------------------------------------------------------------------*/

#include "ecs_descr.h"
#include "ecs_descr_chaine.h"
#include "ecs_maillage.h"
#include "ecs_maillage_priv.h"
#include "ecs_table_att.h"


/*----------------------------------------------------------------------------
 * Fichiers `include' visibles du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_maillage_pre.h"


/*----------------------------------------------------------------------------
 * Fichier `include' du  paquetage courant associé au fichier courant
 *----------------------------------------------------------------------------*/

#include "ecs_pre_cgns.h"


/*----------------------------------------------------------------------------
 * Fichiers `include' privés du paquetage courant
 *----------------------------------------------------------------------------*/


/*============================================================================
 * Définitions de paramètres et macros
 *============================================================================*/

#define ECS_CGNS_TAILLE_NOM    32      /* Longueur nom max (la documentation
                                          ne précise rien, mais les exemples
                                          CGNS et le fichier ADF.h prévoient
                                          des chaînes de 32 caractères) */

#define ECS_CGNS_SSELT_NBR_MAX_SOM   4  /* Nombre maximal de sommets
                                           définissant un sous-élément */


/* Compatibilité avec diverses versions de CGNS */

#if defined(CGNS_SCOPE_ENUMS)
#define CS_CG_ENUM(e) CG_ ## e
#else
#define CS_CG_ENUM(e) e
#endif

#if CGNS_VERSION < 3100
#define cgsize_t int
#endif

/*============================================================================
 * Définitions de structures locales
 *============================================================================*/

/* Définition d'une base CGNS */

typedef struct {

  char                * nom_fic;        /* Nom du fichier CGNS                */
  int                   num_fic;        /* Identificateur de fichier CGNS     */

  int                   num_base;       /* Identificateur de base CGNS        */

  int                   dim_espace;     /* Dimension de l'espace associé      */
  int                   dim_entite;     /* Dimension entités maillage         */

} ecs_loc_cgns_base_t;


/* Définition de conditions aux limites */

typedef struct {
  int              ind_nom;                         /* Numéro de nom associé  */
  int              num_zone;                        /* Num zone associée      */
  CS_CG_ENUM(GridLocation_t)   support;             /* Support associé        */
  CS_CG_ENUM(PointSetType_t)   ptset_type;          /* Type de liste          */
  int              npnts;                           /* Nombre de points       */
  cgsize_t        *pnts;                            /* Liste de points        */
} ecs_loc_cgns_boco_t;


/* Définition d'une section (pour zones non structurées) */

typedef struct {

  char           nom[ECS_CGNS_TAILLE_NOM + 1];      /* Nom de la section      */
  CS_CG_ENUM(ElementType_t)  type;                  /* Type élément           */
  cgsize_t       num_elt_deb;                       /* Numéro premier élement */
  cgsize_t       num_elt_fin;                       /* Numéro dernier élement */
  int            nbr_brd;                           /* Nbr. éléments de bord  */
  int            parent;                            /* 0 si pas de parents,
                                                       1 sinon                */
  cgsize_t      *offsets;                           /* Offsets                */
  cgsize_t      *elems;                             /* Connect. temporaire    */
} ecs_loc_cgns_section_t;


/* Définition d'une zone */

typedef struct {

  char            nom[ECS_CGNS_TAILLE_NOM + 1];   /* Nom de la zone          */
  CS_CG_ENUM(ZoneType_t)      type;               /* Type de la zone         */
  int             nbr_som;                        /* Nombre de sommets       */
  int             nbr_cel;                        /* Nombre de cellules      */
  cgsize_t        taille[3*3];                    /* Indices maximaux pour
                                                     les sommets, cellules,
                                                     et cellules de bord     */
  int             nbr_sections;                   /* Nombre de sections      */
  cgsize_t        num_som_deb;                    /* Numéro premier sommet   */
  cgsize_t        num_elt_deb;                    /* Numéro premier élement  */
  cgsize_t        num_elt_fin;                    /* Numéro dernier élement  */
  CS_CG_ENUM(AngleUnits_t)    angle;              /* Unités angles           */
  ecs_loc_cgns_section_t  *tab_sections;          /* Descriptions sections   */
  cgsize_t        renum_size;                     /* Renumbering array size  */
  cgsize_t       *renum;                          /* Numéro CGNS -> ECS      */

} ecs_loc_cgns_zone_t;


/* Map element type from CGNS to Preprocessor */

typedef struct {

  ecs_elt_typ_t               ecs_type;    /* ECS element type               */
  ecs_int_t                   nbr_som;     /* Number of matching vertices    */
  ecs_int_t                   num_som[8];  /* ECS vertices */

} ecs_cgns_elt_t;


/*============================================================================
 * Définitions de variables globales
 *============================================================================*/


/*============================================================================
 * Fonctions privées
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Création d'une structure `ecs_loc_cgns_base_t' et ouverture d'un
 * fichier CGNS en lecture
 *----------------------------------------------------------------------------*/

static ecs_loc_cgns_base_t *
ecs_loc_pre_cgns__cree(const char  *nom_fichier,
                       int          num_maillage)
{

  ecs_loc_cgns_base_t  * base;

  char  nom_tmp[ECS_CGNS_TAILLE_NOM + 1];
  int   nbases, cell_dim, phys_dim;
  int   ind;
  int   ret;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Base et fichier associés */

  ECS_MALLOC(base, 1, ecs_loc_cgns_base_t);

  ECS_MALLOC(base->nom_fic, strlen(nom_fichier) + 1, char);
  strcpy(base->nom_fic, nom_fichier);

  ret = cg_open(base->nom_fic, CG_MODE_READ, &(base->num_fic));

  if (ret < 0)
    ecs_error(__FILE__, __LINE__, 0,
              _("CGNS: error opening file \"%s\":\n%s"),
              nom_fichier, cg_get_error());


  /* Vérification du nombre de bases (affichage nom si plusieurs) */

  if (cg_nbases(base->num_fic, &nbases) != 0)
    ecs_error
      (__FILE__, __LINE__, 0,
       _("CGNS error:\n%s\n\n"
         "This may be due to a CGNS library supporting only ADF (legacy)\n"
         "or HDF5 (future default) and the file not matching. In this case,\n"
         "running the adf2hdf of hdf2adf converters may solve the issue."),
       cg_get_error());

  if (nbases > 1) {
    printf(_("\n  The file contains multiple bases:\n"));
    for (ind = 0; ind < nbases; ind++) {
      if (   cg_base_read(base->num_fic, ind + 1, nom_tmp, &cell_dim, &phys_dim)
          != 0)
        ecs_error(__FILE__, __LINE__, 0,
                  _("CGNS error:\n%s"), cg_get_error());
      printf(_("    Base %2d: %s\n"), ind + 1, nom_tmp);
    }
    if (num_maillage == 0)
      printf(_("\n  No base was requested; "
               "the first one is used\n\n"));
    else if (num_maillage > 0)
      printf(_("\n  Base number %d was requested\n\n"), num_maillage);
  }
  else if (nbases == 0)
    ecs_error(__FILE__, __LINE__, 0,
              _("no CGNS base in file.\n"));

  assert (num_maillage >= 0);
  if (num_maillage > nbases)
    ecs_error(__FILE__, __LINE__, 0,
              _("The mesh number requested (%d) is higher than\n"
                "the number of bases defined (%d) in file\n%s.\n"),
              num_maillage, nbases, nom_fichier);
  else
    base->num_base = ECS_MAX(1, num_maillage);


  /* Fin d'initialisation et renvoi de la structure */

  base->dim_entite = 0;
  base->dim_espace = 0;

  return base;
}

/*----------------------------------------------------------------------------
 * Fermeture d'un fichier CGNS en lecture et destruction de la structure
 * `ecs_loc_cgns_base_t' associée
 *----------------------------------------------------------------------------*/

static ecs_loc_cgns_base_t *
ecs_loc_pre_cgns__detruit(ecs_loc_cgns_base_t  *base)
{

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(base != NULL);

  if (cg_close(base->num_fic) != 0)
    ecs_error(__FILE__, __LINE__, 0,
              _("CGNS: error closing file \"%s\":\n%s"),
              base->nom_fic, cg_get_error());

  ECS_FREE(base->nom_fic);

  ECS_FREE(base);

  return base;
}

/*----------------------------------------------------------------------------
 * Affichage du titre associé à une base CGNS
 *----------------------------------------------------------------------------*/

static void
ecs_loc_pre_cgns__aff_titre_cas(int          num_fic,
                                int          num_base,
                                const char  *nom_rub)
{
  char  nom_tmp[ECS_CGNS_TAILLE_NOM + 1];
  int   n_desc, ind_desc;
  int   ind_cur, ind_fin, ind_deb, l_dec;

  const int l_ligne = 72;

  char *text = NULL;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (   cg_goto(num_fic, num_base, "end") == 0
      && cg_ndescriptors(&n_desc) == 0) {

    if (n_desc > 0) {
      for (ind_desc = 0; ind_desc < n_desc; ind_desc++) {

        /*
          L'appel cg_descriptor_read alloue text;
          on devra désallouer après utilisation
        */

        if (   cg_descriptor_read(ind_desc + 1, nom_tmp, &text) == 0
            && strcmp(nom_tmp, "CaseTitle") == 0) {

          printf("%s", nom_rub);

          l_dec = strlen(nom_rub);
          ind_fin = strlen(text);
          ind_deb   = 0;
          ind_cur   = l_ligne - l_dec;

          /*
            Découpage de l'affichage pour placer des retours à la
            ligne sur titres longs
          */

          while (ind_cur < ind_fin) {

            while (   ind_cur > ind_deb
                   && text[ind_cur] != ' '
                   && text[ind_cur] != '\n'
                   && text[ind_cur] != '\t')
              ind_cur--;

            if (ind_cur > ind_deb)
              text[ind_cur] = '\0';
            else
              break;

            printf("%s\n%*s", text + ind_deb, l_dec, " ");

            ind_deb = ind_cur + 1;
            while (   ind_deb < ind_fin
                   && text[ind_deb] == ' '
                   && text[ind_deb] == '\n'
                   && text[ind_deb] == '\t')
              ind_deb++;

            ind_cur = ind_deb + l_ligne - l_dec;

          }

          if (ind_deb < ind_fin)
            printf("%s\n", text + ind_deb);

        }
        free(text);
      }

    }
  }

}

/*----------------------------------------------------------------------------
 * Lecture d'informations sur les zones
 *----------------------------------------------------------------------------*/

static ecs_loc_cgns_zone_t *
ecs_loc_pre_cgns__lit_zones(const ecs_loc_cgns_base_t  *base_maillage,
                            int                         nzones)
{
  /* Déclarations des variables de stockage        */
  /* avant transfert dans la structure du maillage */
  /*-----------------------------------------------*/

  char          nom_tmp [ECS_CGNS_TAILLE_NOM + 1];

  const char   *nom_coo_type[] = {N_("cartesian"),
                                  N_("cylindrical"),
                                  N_("spherical")};

  ecs_int_t     coo_type;

  ecs_int_t     nbr_cel_loc;
  ecs_int_t     nbr_som_loc;
  ecs_int_t     nbr_som_tot;

  bool          connect_1to1;
  bool          connect_abutt;
  bool          connect_overset;

  ecs_loc_cgns_zone_t     *ptr_zone;
  ecs_loc_cgns_section_t  *ptr_section;

  ecs_loc_cgns_zone_t     *tab_zone;

  /* Déclarations des variables pour CGNS */
  /*-------------------------------------*/

  int         ind_conn;
  int         ind_section;
  int         ind_zone;
  int         nconns;
  int         num_fic;
  int         num_base;
  int         num_section;
  int         num_zone;
  int         ngrids;

  CS_CG_ENUM(MassUnits_t)         mass;
  CS_CG_ENUM(LengthUnits_t)       length;
  CS_CG_ENUM(TimeUnits_t)         time;
  CS_CG_ENUM(TemperatureUnits_t)  temperature;
  CS_CG_ENUM(AngleUnits_t)        angle;

  CS_CG_ENUM(DataType_t)  type_coord_lu;

  int         ret = 0;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Initialisations */

  num_fic  = base_maillage->num_fic;
  num_base = base_maillage->num_base;

  ECS_MALLOC(tab_zone, nzones, ecs_loc_cgns_zone_t);

  nbr_som_tot = 0;

  connect_1to1    = false;
  connect_abutt   = false;
  connect_overset = false;

  /* Boucle sur les zones */
  /*----------------------*/

  for (ind_zone = 0; ind_zone < nzones; ind_zone++) {

    num_zone = ind_zone + 1;
    ptr_zone = tab_zone + ind_zone;

    ptr_zone->renum_size = -1;
    ptr_zone->num_elt_deb = INT_MAX;
    ptr_zone->num_elt_fin = INT_MIN;
    ptr_zone->renum = NULL;

    /* Informations principales sur la zone */

    if (   cg_zone_read(num_fic, num_base, num_zone, ptr_zone->nom,
                        ptr_zone->taille) != CG_OK
        || cg_zone_type(num_fic, num_base, num_zone,
                        &(ptr_zone->type)) != CG_OK)
      ecs_error(__FILE__, __LINE__, 0,
                _("CGNS error:\n%s"), cg_get_error());

    /* Informations sur le type de coordonnées */

    coo_type = 0;

    ret = cg_coord_info(num_fic, num_base, num_zone, 1, &type_coord_lu,
                        nom_tmp);

    if (ret == CG_OK) {

      if (strcmp("CoordinateR", nom_tmp) == 0)
        coo_type = 1;

      if (base_maillage->dim_espace == 3) {
        ret = cg_coord_info(num_fic, num_base, num_zone, 3, &type_coord_lu,
                            nom_tmp);
        if (   ret == CG_OK
            && (coo_type > 0 && strcmp("CoordinatePhi", nom_tmp) == 0))
          coo_type = 2;
      }

    }
    ret = cg_goto(num_fic, num_base, "Zone_t", num_zone,
                  "GridCoordinates_t", 1, "end");

    if (ret != CG_OK)
      ecs_error(__FILE__, __LINE__, 0,
                _("CGNS error:\n%s"), cg_get_error());

    /* Type d'unité (information optionnelle) */

    length = CS_CG_ENUM(LengthUnitsNull);
    ret = cg_goto(num_fic, num_base, "Zone_t", num_zone,
                  "GridCoordinates_t", 1, "end");
    if (ret == CG_OK)
      ret = cg_units_read(&mass, &length, &time, &temperature, &angle);

    if (ret == CG_OK)
      ptr_zone->angle = angle;

    /* Affichage et stockage des dimensions */

    nbr_som_loc = 0;
    nbr_cel_loc = 0;

    printf(_("\n    Zone %d: \"%s\"; type \"%s\"\n"),
           num_zone, ptr_zone->nom, ZoneTypeName[ptr_zone->type]);

    if (   ptr_zone->type == CS_CG_ENUM(Structured)
        && base_maillage->dim_entite > 1) {
      if (base_maillage->dim_entite == 2) {
        nbr_som_loc = ptr_zone->taille[0] * ptr_zone->taille[1];
        nbr_cel_loc = ptr_zone->taille[2] * ptr_zone->taille[3];

        printf(_("\n      %lu (%lu x %lu) vertices; %lu (%lu x %lu) cells\n"),
               (unsigned long)nbr_som_loc,
               (unsigned long)ptr_zone->taille[0],
               (unsigned long)ptr_zone->taille[1],
               (unsigned long)nbr_cel_loc,
               (unsigned long)ptr_zone->taille[2],
               (unsigned long)ptr_zone->taille[3]);
      }
      else if (base_maillage->dim_entite == 3) {
        nbr_som_loc =   ptr_zone->taille[0] * ptr_zone->taille[1]
                      * ptr_zone->taille[2];
        nbr_cel_loc =   ptr_zone->taille[3] * ptr_zone->taille[4]
                      * ptr_zone->taille[5];

        printf(_("\n      %lu (%lu x %lu x %lu) vertices;"
                 " %lu (%lu x %lu x %lu) cells\n"),
               (unsigned long)nbr_som_loc,
               (unsigned long)ptr_zone->taille[0],
               (unsigned long)ptr_zone->taille[1],
               (unsigned long)ptr_zone->taille[2],
               (unsigned long)nbr_cel_loc,
               (unsigned long)ptr_zone->taille[3],
               (unsigned long)ptr_zone->taille[4],
               (unsigned long)ptr_zone->taille[5]);
      }
    }
    else {
      nbr_som_loc = ptr_zone->taille[0];
      nbr_cel_loc = ptr_zone->taille[1];

      printf(_("\n      %lu vertices; %lu cells\n"),
             (unsigned long)nbr_som_loc, (unsigned long)nbr_cel_loc);
    }

    ptr_zone->nbr_som = nbr_som_loc;
    ptr_zone->nbr_cel = nbr_cel_loc;

    ptr_zone->num_som_deb  = nbr_som_tot + 1;
    nbr_som_tot           += nbr_som_loc;

    /* Informations complémentaires (non indispensables) */

    printf(_("      (%s coordinates, \"%s\" precision, unit \"%s\")\n"),
           nom_coo_type[coo_type],
           DataTypeName[type_coord_lu], LengthUnitsName[length]);

    if (cg_ngrids(num_fic, num_base, num_zone, &ngrids) == CG_OK) {
      if (ngrids > 1)
        printf(_("      %d time-varying coordinates "
                 "(the first are used)\n"), ngrids);
    }

    /* Informations sur les connectivités multizones */

    ret = cg_n1to1(num_fic, num_base, num_zone, &nconns);

    if (ret == CG_OK) {
      if (nconns > 0)
        connect_1to1 = true;
    }

    ret = cg_nconns(num_fic, num_base, num_zone, &nconns);

    if (ret == CG_OK) {

      char nom_tmp_aux [ECS_CGNS_TAILLE_NOM + 1];

      CS_CG_ENUM(GridLocation_t)         location;
      CS_CG_ENUM(GridConnectivityType_t) connect_type;
      CS_CG_ENUM(PointSetType_t)         ptset_type;
      cgsize_t                           npnts;
      CS_CG_ENUM(ZoneType_t)             donor_zonetype;
      CS_CG_ENUM(PointSetType_t)         donor_ptset_type;
      CS_CG_ENUM(DataType_t)             donor_datatype;
      cgsize_t                           ndata_donor;

      for (ind_conn = 0; ind_conn < nconns; ind_conn++) {

        ret = cg_conn_info(num_fic, num_base, num_zone, ind_conn + 1,
                           nom_tmp, &location, &connect_type, &ptset_type,
                           &npnts, nom_tmp_aux, &donor_zonetype,
                           &donor_ptset_type, &donor_datatype,
                           &ndata_donor);

        if (ret == CG_OK) {

          switch (connect_type) {
          case CS_CG_ENUM(Overset):
            connect_overset = true;
            break;
          case CS_CG_ENUM(Abutting):
            connect_abutt = true;
            break;
          case CS_CG_ENUM(Abutting1to1):
            connect_1to1 = true;
            break;
          default:
            break;
          }

        }
      }
    }

    /* Informations sur le nombre de sections (non structuré) */

    ptr_zone->nbr_sections = 0;
    ptr_zone->tab_sections = NULL;

    if (ptr_zone->type == CS_CG_ENUM(Structured)) {

      ptr_zone->num_elt_deb = 1;
      ptr_zone->num_elt_fin = ptr_zone->nbr_cel;

    }
    else {

      ret = cg_nsections(num_fic, num_base, num_zone,
                         &(ptr_zone->nbr_sections));
      if (ret == CG_OK) {

        ECS_MALLOC(ptr_zone->tab_sections, ptr_zone->nbr_sections,
                   ecs_loc_cgns_section_t);

        for (ind_section = 0;
             ind_section < ptr_zone->nbr_sections;
             ind_section++) {

          num_section = ind_section + 1;
          ptr_section = ptr_zone->tab_sections + ind_section;

          if (cg_section_read(num_fic,
                              num_base,
                              num_zone,
                              num_section,
                              ptr_section->nom,
                              &(ptr_section->type),
                              &(ptr_section->num_elt_deb),
                              &(ptr_section->num_elt_fin),
                              &(ptr_section->nbr_brd),
                              &(ptr_section->parent)) != CG_OK)
            ecs_error(__FILE__, __LINE__, 0,
                      _("CGNS error:\n%s"), cg_get_error());

          ptr_section->offsets = NULL;
          ptr_section->elems = NULL;

          printf(_("      Section %2d: \"%s\";\n"
                   "                   (indices %lu to %lu, type \"%s\")\n"),
                 num_section, ptr_section->nom,
                 (unsigned long)ptr_section->num_elt_deb,
                 (unsigned long)ptr_section->num_elt_fin,
                 ElementTypeName[ptr_section->type]);

          ptr_zone->num_elt_deb = ECS_MIN(ptr_section->num_elt_deb,
                                          ptr_zone->num_elt_deb);
          ptr_zone->num_elt_fin = ECS_MAX(ptr_section->num_elt_fin,
                                          ptr_zone->num_elt_fin);

        }

      }

    }

  }

  printf("\n");

  /* Informations sur les connectivités */

  if (connect_1to1 == true) {
    ecs_warn();
    printf(_("The CGNS mesh read contains multizone (\"one to one\")\n"
             "vertex equivalences which are not automatically handled\n"
             "by the Preprocessor.\n"
             "-> Use an appropriate joining option\n"));
  }

  if (connect_abutt == true) {
    ecs_warn();
    printf(_("The CGNS mesh read contains non-conforming (\"abutting\")\n"
             "connectivity information which is not automatically handled\n"
             "by the Preprocessor.\n"
             "-> Use an appropriate joining option\n"));
  }

  if (connect_overset == true) {
    ecs_warn();
    printf(_("The CGNS mesh read contains (\"overset\") connectivity\n"
             "information which is not handled by the Preprocessor.\n"));
  }

  /* Renvoi du tableau d'information sur les zones */

  return tab_zone;
}

/*----------------------------------------------------------------------------
 * Lecture des conditions aux limites
 *----------------------------------------------------------------------------*/

static void
ecs_loc_pre_cgns__lit_boco(const ecs_loc_cgns_base_t    *base_maillage,
                           bool                          ignore_vertex_bocos,
                           const int                     nzones,
                           ecs_loc_cgns_zone_t          *tab_zone,
                           int                          *nbr_nom_boco,
                           int                          *nbr_boco_tot,
                           char                       ***tab_nom_boco,
                           ecs_loc_cgns_boco_t         **tab_boco)
{

  char    nom_tmp [ECS_CGNS_TAILLE_NOM + 1];
  char    nom_fam [ECS_CGNS_TAILLE_NOM + 1];

  int    *nbr_boco_loc;

  /* Déclarations des variables de stockage        */
  /* avant transfert dans la structure du maillage */
  /*-----------------------------------------------*/

  ecs_loc_cgns_zone_t     *ptr_zone;

  ecs_loc_cgns_boco_t     *tab_boco_loc = NULL;
  char                   **nom_boco_loc = NULL;

  /* Déclarations des variables pour CGNS */
  /*-------------------------------------*/

  int         ind_boco;
  int         ind_boco_glob;
  int         ind_nom;
  int         ind_zone;
  int         num_boco;
  int         num_fic;
  int         num_base;
  int         num_zone;

  int         nbocos;
  cgsize_t    npnts;

  CS_CG_ENUM(BCType_t)        bocotype;
  CS_CG_ENUM(GridLocation_t)  GridLocation;
  CS_CG_ENUM(PointSetType_t)  ptset_type;
  int                         NormalIndex[3];
  cgsize_t                    NormalListFlag;
  CS_CG_ENUM(DataType_t)      NormalDataType;
  int                         ndataset;

  void        *normales;

  int         ret = CG_OK;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Initialisations */

  num_fic  = base_maillage->num_fic;
  num_base = base_maillage->num_base;

  ECS_MALLOC(nbr_boco_loc, nzones, int);

  *nbr_nom_boco = 0;
  *tab_nom_boco = NULL;

  *tab_boco = NULL;

  /* Comptage */

  *nbr_boco_tot = 0;

  for (ind_zone = 0; ind_zone < nzones; ind_zone++) {

    num_zone = ind_zone + 1;
    ptr_zone = tab_zone + ind_zone;

    ret = cg_nbocos(num_fic, num_base, num_zone, &nbocos);

    if (ret == CG_OK) {
      nbr_boco_loc[ind_zone] = nbocos;
      *nbr_boco_tot += nbocos;
    }
    else
      ecs_error(__FILE__, __LINE__, 0,
                _("CGNS error:\n%s"), cg_get_error());

  }

  if (*nbr_boco_tot > 0) {

    printf(_("  CGNS boundary condition information:\n"));

    ECS_MALLOC(tab_boco_loc, *nbr_boco_tot, ecs_loc_cgns_boco_t);
    ECS_MALLOC(nom_boco_loc, *nbr_boco_tot, char *);

    for (ind_nom = 0; ind_nom < *nbr_boco_tot; ind_nom++)
      nom_boco_loc[ind_nom] = NULL;

    *tab_boco = tab_boco_loc;

  }

  /*------------------------------------*/
  /* Lecture des conditions aux limites */
  /*------------------------------------*/

  ind_boco_glob = 0;

  for (ind_zone = 0; ind_zone < nzones; ind_zone++) {

    num_zone = ind_zone + 1;
    ptr_zone = tab_zone + ind_zone;

    if (nbr_boco_loc[ind_zone] > 0) {

      printf(_("\n    Zone %d\n"), num_zone);

      for (ind_boco = 0; ind_boco < nbr_boco_loc[ind_zone]; ind_boco++) {

        num_boco = ind_boco + 1;

        /* Informations principales */

        ret = cg_boco_info(num_fic, num_base, num_zone, num_boco, nom_tmp,
                           &bocotype, &ptset_type, &npnts, NormalIndex,
                           &NormalListFlag, &NormalDataType, &ndataset);

        if (ret != CG_OK)
          break;

        nom_tmp [ECS_CGNS_TAILLE_NOM] = '\0';

        /* Informations sur support et impression */

        GridLocation = CS_CG_ENUM(Vertex);

        ret = cg_goto(num_fic, num_base, "Zone_t", num_zone,
                      "ZoneBC_t", 1, "BC_t", num_boco, "end");

        if (ret == CG_OK) {

          ret = cg_famname_read(nom_fam);

          if (ret == CG_OK) {
            nom_fam[ECS_CGNS_TAILLE_NOM] = '\0';
            strcpy(nom_tmp, nom_fam);
          }

          if (ptr_zone->type != CS_CG_ENUM(Structured)) {

            ret = cg_gridlocation_read(&GridLocation);

            if (ret != CG_OK)
              GridLocation = CS_CG_ENUM(Vertex);

          }

        }

        if (ptr_zone->type == CS_CG_ENUM(Structured)) {

          printf(_("      BC %2d: \"%s\" (\"%s\")\n"),
                 num_boco, nom_tmp, BCTypeName[bocotype]);

        }
        else {

          if (   ptset_type == CS_CG_ENUM(ElementList)
              || ptset_type == CS_CG_ENUM(ElementRange))
            GridLocation = CS_CG_ENUM(FaceCenter);

          printf(_("      BC %2d: \"%s\" (\"%s\" on \"%s\")\n"),
                 num_boco, nom_tmp, BCTypeName[bocotype],
                 GridLocationName[GridLocation]);

          if (GridLocation != CS_CG_ENUM(Vertex)) {
            if (ptr_zone->renum_size < 0)
              ptr_zone->renum_size = 0;
          }

        }

        if (ignore_vertex_bocos && GridLocation == CS_CG_ENUM(Vertex))
          continue;

        /* stockage */

        for (ind_nom = 0; nom_boco_loc[ind_nom] != NULL; ind_nom++) {
          if (strcmp(nom_boco_loc[ind_nom], nom_tmp) == 0)
            break;
        }
        if (nom_boco_loc[ind_nom] == NULL) {
          ECS_MALLOC(nom_boco_loc[ind_nom], strlen(nom_tmp) + 1, char);
          strcpy(nom_boco_loc[ind_nom], nom_tmp);
        }

        (tab_boco_loc[ind_boco_glob]).ind_nom     = ind_nom;
        (tab_boco_loc[ind_boco_glob]).num_zone    = num_zone;
        (tab_boco_loc[ind_boco_glob]).support     = GridLocation;
        (tab_boco_loc[ind_boco_glob]).ptset_type  = ptset_type;
        (tab_boco_loc[ind_boco_glob]).npnts       = npnts;

        /* Lecture listes */

        if (   ptset_type == CS_CG_ENUM(PointRange)
            || ptset_type == CS_CG_ENUM(ElementRange))
          ECS_MALLOC((tab_boco_loc[ind_boco_glob]).pnts, npnts * 3, cgsize_t);
        else
          ECS_MALLOC((tab_boco_loc[ind_boco_glob]).pnts, npnts, cgsize_t);

        if (NormalListFlag > 0) {
          if (NormalDataType == CS_CG_ENUM(RealSingle))
            ECS_MALLOC(normales, NormalListFlag, float);
          else
            ECS_MALLOC(normales, NormalListFlag, double);
        }
        else
          normales = NULL;

        ret = cg_boco_read(num_fic, num_base, num_zone, num_boco,
                           (tab_boco_loc[ind_boco_glob]).pnts,
                           normales);

        if (ret != CG_OK)
          break;

        if (normales != NULL)
          ECS_FREE(normales);

        /* Incrémentation compteur */

        ind_boco_glob += 1;

      }

      if (ret != CG_OK)
        ecs_error(__FILE__, __LINE__, 0,
                  _("CGNS error:\n%s"), cg_get_error());

    }

  }

  if (*nbr_boco_tot > 0)
    printf("\n");

  if (ind_boco_glob < *nbr_boco_tot)
    *nbr_boco_tot = ind_boco_glob;

  if (*nbr_boco_tot > 0) {

    for (ind_nom = 0;
         ind_nom < *nbr_boco_tot && nom_boco_loc[ind_nom] != NULL;
         ind_nom++);

    ECS_REALLOC(nom_boco_loc, ind_nom, char *);

    *nbr_nom_boco = ind_nom;
    *tab_nom_boco = nom_boco_loc;

  }

  ECS_FREE(nbr_boco_loc);
}

/*----------------------------------------------------------------------------
 * Lecture des sommets
 *----------------------------------------------------------------------------*/

static void
ecs_loc_pre_cgns__lit_som(ecs_maillage_t             *maillage,
                          const ecs_loc_cgns_base_t  *base_maillage,
                          int                         nzones,
                          ecs_loc_cgns_zone_t        *tab_zone)
{
  /* Déclarations des variables de stockage        */
  /* avant transfert dans la structure du maillage */
  /*-----------------------------------------------*/

  ecs_int_t     coo_type;
  ecs_int_t     ind_som;
  ecs_int_t     nbr_som;
  ecs_coord_t  *som_val_tmp;
  ecs_coord_t  *som_val_coord;                    /* Coordonnées des sommets  */

  double        cnv_angle;
  double        r_tmp;
  double        t_tmp;
  double        p_tmp;

  ecs_loc_cgns_zone_t  *ptr_zone;

  /* Déclarations des variables pour CGNS */
  /*-------------------------------------*/

  int         ind_zone;
  int         num_fic;
  int         num_base;
  int         num_zone;
  int         phys_dim;
  cgsize_t    irmin[3];
  cgsize_t    irmax[3];

  int        *ind_som_deb;

  CS_CG_ENUM(DataType_t)  type_coord;

  int         ret = 0;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Initialisations */

  num_fic  = base_maillage->num_fic;
  num_base = base_maillage->num_base;

  phys_dim = base_maillage->dim_espace;

  if (sizeof(ecs_coord_t) == sizeof(float))
    type_coord = CS_CG_ENUM(RealSingle);
  else if (sizeof(ecs_coord_t) == sizeof(double))
    type_coord = CS_CG_ENUM(RealDouble);
  else
    assert (   (sizeof(ecs_coord_t) == sizeof(float))
            || (sizeof(ecs_coord_t) == sizeof(double)));

  ECS_MALLOC(ind_som_deb, nzones + 1, int);

  ind_som_deb[0] = 0;

  /* Dimensionnement pour lectures */
  /*-------------------------------*/

  for (ind_zone = 0; ind_zone < nzones; ind_zone++) {

    num_zone = ind_zone + 1;
    ptr_zone = tab_zone + ind_zone;

    ind_som_deb[ind_zone + 1] = ind_som_deb[ind_zone] + ptr_zone->nbr_som;

  }

  ECS_MALLOC(som_val_tmp,   ind_som_deb[nzones] * phys_dim,  ecs_coord_t);
  ECS_MALLOC(som_val_coord, ind_som_deb[nzones] * 3, ecs_coord_t);

  nbr_som = ind_som_deb[nzones];

  /* Lecture effective des coordonnées et éléments*/
  /*----------------------------------------------*/

  for (ind_zone = 0; ind_zone < nzones; ind_zone++) {

    num_zone = ind_zone + 1;
    ptr_zone = tab_zone + ind_zone;

    /*
      Lecture des coordonnées (CGNS convertit automatiquement du type
      simple précision ou double précision au type demandé si nécessaire).
    */

    /* Cas structuré */

    if (ptr_zone->type == CS_CG_ENUM(Structured)) {
      irmin[0] = 1;
      irmin[1] = 1;
      irmax[0] = ptr_zone->taille[0];
      irmax[1] = ptr_zone->taille[1];
      if (phys_dim == 3) {
        irmin[2] = 1;
        irmax[2] = ptr_zone->taille[2];
      }
    }

    /* Cas non structuré */

    else if (ptr_zone->type == CS_CG_ENUM(Unstructured)) {
      irmin[0] = 1;
      irmin[1] = 1;
      irmax[0] = ptr_zone->taille[0];
      irmax[1] = 1;
      irmin[2] = 1;
      irmax[2] = 1;
    }

    /* Première coordonnée (X ou R) */

    coo_type = 0;

    ret = cg_coord_read(num_fic, num_base, num_zone, "CoordinateX",
                        type_coord, irmin, irmax,
                        som_val_tmp + ind_som_deb[ind_zone]);
    if (ret != 0) {
      ret = cg_coord_read(num_fic, num_base, num_zone, "CoordinateR",
                          type_coord, irmin, irmax,
                          som_val_tmp + ind_som_deb[ind_zone] + nbr_som);
      if (ret == 0)
        coo_type = 1;
    }

    /* Seconde coordonnée (Y ou Theta) */

    if (ret == 0 && coo_type == 0)
      ret = cg_coord_read(num_fic, num_base, num_zone, "CoordinateY",
                          type_coord, irmin, irmax,
                          som_val_tmp + ind_som_deb[ind_zone] + nbr_som);
    else if (ret == 0 && coo_type == 1)
        ret = cg_coord_read(num_fic, num_base, num_zone, "CoordinateTheta",
                          type_coord, irmin, irmax,
                          som_val_tmp + ind_som_deb[ind_zone] + nbr_som);

    /* Troisième coordonnée (Z ou Phi) */

    if (ret == 0 && phys_dim == 3) {
      ret = cg_coord_read(num_fic, num_base, num_zone, "CoordinateZ",
                          type_coord, irmin, irmax,
                          som_val_tmp + ind_som_deb[ind_zone] + (nbr_som*2));
      if (ret != 0 && coo_type == 1) {
        ret = cg_coord_read(num_fic, num_base, num_zone, "CoordinatePhi",
                          type_coord, irmin, irmax,
                          som_val_tmp + ind_som_deb[ind_zone] + nbr_som);
        if (ret == 0)
          coo_type = 2;
      }
    }

    if (ret != 0)
      ecs_error(__FILE__, __LINE__, 0,
                _("CGNS error:\n%s"), cg_get_error());

    /* Conversion en coordonnées Cartésiennes si nécessaire */

    if (coo_type > 0) {
      if (ptr_zone->angle == CS_CG_ENUM(Degree))
        cnv_angle = 4 * atan(1) / 180.0;
      else if (ptr_zone->angle == CS_CG_ENUM(Radian))
        cnv_angle = 1.0;
      else {
        cnv_angle = 1.0;
        ecs_warn();
        printf(_("Cylindrical or spherical coordinates with unknown\n"
                 "or undefined angle unit "
                 "(-> radians are considered)\n"));
      }
    }

    if (coo_type == 1 && phys_dim > 1) {
      for (ind_som = 0; ind_som < nbr_som; ind_som++) {
        r_tmp = som_val_tmp[ind_som           ];
        t_tmp = som_val_tmp[ind_som + nbr_som ] * cnv_angle;
        som_val_tmp[ind_som          ] = r_tmp * cos(t_tmp);
        som_val_tmp[ind_som + nbr_som] = r_tmp * sin(t_tmp);
      }
    }
    else if (coo_type == 2 && phys_dim > 2) {
      for (ind_som = 0; ind_som < nbr_som; ind_som++) {
        r_tmp = som_val_tmp[ind_som                ];
        t_tmp = som_val_tmp[ind_som +      nbr_som ] * cnv_angle;
        p_tmp = som_val_tmp[ind_som + (2 * nbr_som)] * cnv_angle;
        som_val_tmp[ind_som                ] = r_tmp * sin(t_tmp) * cos(p_tmp);
        som_val_tmp[ind_som +      nbr_som ] = r_tmp * sin(t_tmp) * sin(p_tmp);
        som_val_tmp[ind_som + (2 * nbr_som)] = r_tmp * cos(t_tmp);
      }
    }

  }

  /* Conversion au format local */

  if (phys_dim == 1) {
    for (ind_som = 0; ind_som < nbr_som; ind_som++) {
      som_val_coord[ind_som * 3    ] = som_val_tmp[ind_som];
      som_val_coord[ind_som * 3 + 1] = 0.0;
      som_val_coord[ind_som * 3 + 2] = 0.0;
    }
  }
  else if (phys_dim == 2) {
    for (ind_som = 0; ind_som < nbr_som; ind_som++) {
      som_val_coord[ind_som * 3    ] = som_val_tmp[ind_som];
      som_val_coord[ind_som * 3 + 1] = som_val_tmp[ind_som + nbr_som ];
      som_val_coord[ind_som * 3 + 1] = 0.0;
    }
  }
  else if (phys_dim == 3) {
    for (ind_som = 0; ind_som < nbr_som; ind_som++) {
      som_val_coord[ind_som * 3    ] = som_val_tmp[ind_som                ];
      som_val_coord[ind_som * 3 + 1] = som_val_tmp[ind_som +      nbr_som ];
      som_val_coord[ind_som * 3 + 2] = som_val_tmp[ind_som + (2 * nbr_som)];
    }
  }
  else {
    ECS_FREE(som_val_coord);
  }

  ECS_FREE(som_val_tmp);

  /* Transfert des valeurs lues dans la structure d'entité de maillage */
  /*===================================================================*/

  ecs_maillage_pre__cree_som(maillage,
                             nbr_som,
                             som_val_coord);

  /* Libération des tableaux de travail */

  ECS_FREE(som_val_tmp);
  ECS_FREE(ind_som_deb);
}

/*----------------------------------------------------------------------------
 * Fonction qui marque les sommets correspondant à une liste.
 *----------------------------------------------------------------------------*/

static void
ecs_loc_pre_cgns__marque_som(const ecs_loc_cgns_base_t  *base_maillage,
                             int                         num_marque,
                             ecs_loc_cgns_zone_t        *ptr_zone,
                             CS_CG_ENUM(PointSetType_t)  ptset_type,
                             int                         npnts,
                             cgsize_t                   *pnts,
                             ecs_int_t                  *indic_som)
{

  ecs_int_t    ient_max;

  ecs_int_t    ind_som;
  ecs_int_t    ind_som_deb;

  ecs_int_t    ii, jj, kk;
  ecs_int_t    ni, nj;

  int          cel_dim;


  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


  /* Initialisations */

  cel_dim = base_maillage->dim_entite;

  if (cel_dim == 1)
    return;
  else if (cel_dim == 2)
    ient_max = ECS_ENTMAIL_FAC;
  else
    ient_max = ECS_ENTMAIL_CEL;


  /* Marquage des sommets */
  /*----------------------*/

  ind_som_deb = ptr_zone->num_som_deb - 1;


  /* Points définis par leur étendue */

  if (ptset_type == CS_CG_ENUM(PointRange)) {

    /* Cas structuré */

    if (ptr_zone->type == CS_CG_ENUM(Structured)) {

      if (ient_max == ECS_ENTMAIL_FAC) {

        ni = ptr_zone->taille[0];

        for (jj = pnts[1] - 1; jj < pnts[4]; jj++) {
          for (ii = pnts[0] - 1; ii < pnts[3]; ii++) {
            ind_som = ii + jj*ni;
            indic_som[ind_som_deb + ind_som] = num_marque;
          }
        }

      }
      else if (ient_max == ECS_ENTMAIL_CEL) {

        ni = ptr_zone->taille[0];
        nj = ptr_zone->taille[1];

        for (kk = pnts[2] - 1; kk < pnts[5]; kk++) {
          for (jj = pnts[1] - 1; jj < pnts[4]; jj++) {
            for (ii = pnts[0] - 1; ii < pnts[3]; ii++) {
              ind_som = ii + jj*ni + kk*ni*nj;
              indic_som[ind_som_deb + ind_som] = num_marque;
            }
          }
        }
      }

    }

    /* Cas non structuré */

    else {

      for (ind_som = pnts[0]; ind_som < pnts[1]; ind_som++)
        indic_som[ind_som_deb + ind_som - 1] = num_marque;

    }

  }

  /* Points définis par une liste */

  else {

    for (ind_som = 0; ind_som < npnts; ind_som++)
      indic_som[ind_som_deb + pnts[ind_som] - 1] = num_marque;

  }

}


/*----------------------------------------------------------------------------
 * Fonction qui créé des entités de maillage supplémentaires pour porter
 * des conditions aux limites définies sur les sommets.
 *----------------------------------------------------------------------------*/

static void
ecs_loc_pre__cgns__cree_ent_inf_som(const ecs_loc_cgns_base_t  *base_maillage,
                                    size_t                 nbr_elt_ent[],
                                    ecs_size_t            *elt_pos_som_ent[],
                                    ecs_int_t             *elt_val_som_ent[],
                                    int                    nzones,
                                    ecs_loc_cgns_zone_t   *tab_zone,
                                    int                    nbr_boco_tot,
                                    ecs_loc_cgns_boco_t   *tab_boco,
                                    ecs_int_t             *nbr_boco_som,
                                    ecs_int_t            **ind_nom_boco_som,
                                    ecs_int_t            **nbr_sselt_boco_som,
                                    ecs_int_t           ***tab_sselt_boco_som
)
{

  size_t       cpt_sselt;
  size_t       cpt_sselt_boco;
  size_t       cpt_val_sselt;
  size_t       nbr_sselt_ini;
  size_t       nbr_val_sselt_ini;

  int          ind_boco;
  int          ind_boco_sub;
  ecs_int_t    cpt_boco_som;
  size_t       ind_elt;
  size_t       ind_pos_elt;
  size_t       ind_pos_sselt;
  size_t       ind_som;
  size_t       ind_sselt;

  size_t       nbr_elt;
  size_t       nbr_som;
  size_t       nbr_som_elt;
  int          num_boco;
  ecs_int_t    num_def;

  ecs_int_t    typ_elt;

  bool         boco_som;
  bool         bool_cree;

  int          cel_dim;
  int          ind_zone;

  ecs_int_t   *indic_som;
  ecs_int_t   *indic_sselt;
  ecs_int_t   *liste_sselt;
  ecs_int_t   *renum_sselt;

  ecs_size_t  *pos_som_elt;
  ecs_int_t   *val_som_elt;

  ecs_size_t  *pos_som_sselt;
  ecs_int_t   *val_som_sselt;

  ecs_loc_cgns_zone_t   *ptr_zone;
  ecs_loc_cgns_boco_t   *ptr_boco;

  const ecs_elt_typ_t  *typ_geo_base;


  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructionsxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


  /* Initialisations */

  cel_dim = base_maillage->dim_entite;

  if (cel_dim < 3)
    return;

  ind_boco = 0;

  cpt_boco_som = 0;

  *nbr_boco_som = 0;
  *ind_nom_boco_som = NULL;
  *nbr_sselt_boco_som = NULL;
  *tab_sselt_boco_som = NULL;


  /* Détermination de base du type d'élément "classique" */

  assert(cel_dim >=2 && cel_dim <= 3);

  typ_geo_base = ecs_glob_typ_elt[cel_dim - 2];


  /* On vérifie que l'on a quelque chose à faire */
  /*---------------------------------------------*/

  boco_som = false;

  for (ind_boco = 0; ind_boco < nbr_boco_tot; ind_boco++) {

    ptr_boco = tab_boco + ind_boco;

    if (ptr_boco->support == CS_CG_ENUM(Vertex))
      boco_som = true;

  }

  if (boco_som == false)
    return;

  /* Initialisation du marquage des sommets */
  /*----------------------------------------*/

  nbr_som = 0;

  for (ind_zone = 0; ind_zone < nzones; ind_zone++)
    nbr_som += (tab_zone + ind_zone)->nbr_som;

  ECS_MALLOC(indic_som, nbr_som, ecs_int_t);

  for (ind_som = 0; ind_som < nbr_som; ind_som++)
    indic_som[ind_som] = 0;

  for (ind_boco = 0; ind_boco < nbr_boco_tot; ind_boco++) {

    ptr_boco = tab_boco + ind_boco;

    num_boco = ptr_boco->ind_nom + 1;

    if (ptr_boco->support == CS_CG_ENUM(Vertex)) {

      ptr_zone = tab_zone + ptr_boco->num_zone - 1;

      ecs_loc_pre_cgns__marque_som(base_maillage,
                                   num_boco,
                                   ptr_zone,
                                   ptr_boco->ptset_type,
                                   ptr_boco->npnts,
                                   ptr_boco->pnts,
                                   indic_som);

    }

  }

  /* Comptage pour la création des faces */
  /*-------------------------------------*/

  /*
    On créé des faces lorsque tous les sommets
    d'une face d'une cellule appartiennent à une C.L.
    Ceci crée des faces en trop (i.e. toutes les faces intérieures)
    dans le cas d'un maillage 2D extrudé avec une seule couche
    d'éléments, mais évite d'effectuer une boucle sur les éléments
    par C.L. et par zone. On compactera les faces par la suite
    (avant de passer en connectivité descendante)
  */

  nbr_elt     = nbr_elt_ent[ECS_ENTMAIL_CEL];
  pos_som_elt = elt_pos_som_ent[ECS_ENTMAIL_CEL];
  val_som_elt = elt_val_som_ent[ECS_ENTMAIL_CEL];


  /* Boucle sur cellules  */

  cpt_sselt  = 0;
  cpt_val_sselt  = 0;

  for (ind_elt = 0; ind_elt < nbr_elt; ind_elt++) {

    nbr_som_elt = pos_som_elt[ind_elt + 1] - pos_som_elt[ind_elt];

    typ_elt = typ_geo_base[nbr_som_elt];

    ind_pos_elt = pos_som_elt[ind_elt] - 1;

    /* Boucle sur les faces définissant la cellule */

    for (ind_sselt = 0;
         ind_sselt < (size_t)(ecs_fic_elt_typ_liste_c[typ_elt].nbr_sous_elt);
         ind_sselt++) {

      /* Boucle sur les sommets définissant la face */

      ind_som = 0;

      bool_cree = true;

      while(   ind_som < ECS_CGNS_SSELT_NBR_MAX_SOM
            && bool_cree == true
            && (num_def = ecs_fic_elt_typ_liste_c
                            [typ_elt].sous_elt[ind_sselt].som[ind_som]) != 0) {

        if (indic_som[val_som_elt[ind_pos_elt + num_def - 1] - 1] == 0)
          bool_cree = false;

        ind_som++;

      } /* Fin de la boucle sur les sommets définissant la face */

      if (bool_cree == true) {
        cpt_sselt     += 1;
        cpt_val_sselt += ind_som;
      }

    }  /* Fin de la boucle sur les faces d'une cellule */

  } /* Fin de la boucle sur les éléments */


  if (cpt_sselt == 0) {
    ECS_FREE(indic_som);
    return;
  }

  /* Construction des faces */
  /*------------------------*/

  nbr_sselt_ini = nbr_elt_ent[ECS_ENTMAIL_FAC];

  if (nbr_sselt_ini > 0)
    nbr_val_sselt_ini = elt_pos_som_ent[ECS_ENTMAIL_FAC][nbr_sselt_ini] - 1;
  else
    nbr_val_sselt_ini = 0;

  ECS_REALLOC(elt_pos_som_ent[ECS_ENTMAIL_FAC],
              nbr_sselt_ini + cpt_sselt + 1,
              ecs_size_t);
  ECS_REALLOC(elt_val_som_ent[ECS_ENTMAIL_FAC],
              nbr_val_sselt_ini + cpt_val_sselt,
              ecs_int_t);

  nbr_elt_ent[ECS_ENTMAIL_FAC] += cpt_sselt;

  pos_som_sselt = elt_pos_som_ent[ECS_ENTMAIL_FAC];
  val_som_sselt = elt_val_som_ent[ECS_ENTMAIL_FAC];

  if (nbr_sselt_ini == 0)
    pos_som_sselt[0] = 1;

  /* Boucle sur les éléments  */

  cpt_sselt  = nbr_sselt_ini;
  cpt_val_sselt  = nbr_val_sselt_ini;

  for (ind_elt = 0; ind_elt < nbr_elt; ind_elt++) {

    nbr_som_elt = pos_som_elt[ind_elt + 1] - pos_som_elt[ind_elt];

    typ_elt = typ_geo_base[nbr_som_elt];

    ind_pos_elt = pos_som_elt[ind_elt] - 1;

    /* Boucle sur les sous-éléments définissant l'élément */

    for (ind_sselt = 0;
         ind_sselt < (size_t)(ecs_fic_elt_typ_liste_c[typ_elt].nbr_sous_elt);
         ind_sselt++) {

      /* Boucle sur les sommets définissant le sous-élément */

      ind_som = 0;

      bool_cree = true;

      while(   ind_som < ECS_CGNS_SSELT_NBR_MAX_SOM
            && bool_cree == true
            && (num_def = ecs_fic_elt_typ_liste_c
                [typ_elt].sous_elt[ind_sselt].som[ind_som]) != 0) {

        if (indic_som[val_som_elt[ind_pos_elt + num_def - 1] - 1] == 0)
          bool_cree = false;

        ind_som++;

      }


      /* Si l'on crée une face */

      if (bool_cree == true) {

        ind_som = 0;

        while(   ind_som < ECS_CGNS_SSELT_NBR_MAX_SOM
              && bool_cree == true
              && (num_def
                  = ecs_fic_elt_typ_liste_c
                      [typ_elt].sous_elt[ind_sselt].som[ind_som]) != 0) {

          /* Définition de la face en fonction des sommets */

          val_som_sselt[cpt_val_sselt++]
            = val_som_elt[ind_pos_elt + num_def - 1];

          ind_som++;

        }

        /* Position de la face dans sa définition en fonction des sommets */

        pos_som_sselt[cpt_sselt + 1] = pos_som_sselt[cpt_sselt] + ind_som;

        cpt_sselt     += 1;

      } /* Fin création face */

    } /* Fin de la boucle sur les faces d'une cellule */

  } /* Fin de la boucle sur les éléments */

  /* Création listes d'éléments inférieurs et préparation compactage */
  /*-----------------------------------------------------------------*/

  ECS_MALLOC(indic_sselt, nbr_elt_ent[ECS_ENTMAIL_FAC], ecs_int_t);
  ECS_MALLOC(liste_sselt, nbr_elt_ent[ECS_ENTMAIL_FAC], ecs_int_t);

  for (ind_sselt = 0; ind_sselt < nbr_elt_ent[ECS_ENTMAIL_FAC]; ind_sselt++)
    indic_sselt[ind_sselt] = 0;

  cpt_sselt = nbr_sselt_ini;

  for (ind_boco = 0; ind_boco < nbr_boco_tot; ind_boco++) {

    ptr_boco = tab_boco + ind_boco;

    if (ptr_boco->ind_nom > -1 && ptr_boco->support == CS_CG_ENUM(Vertex)) {

      num_boco = ptr_boco->ind_nom + 1;

      ECS_REALLOC(*ind_nom_boco_som,   cpt_boco_som + 1, ecs_int_t);
      ECS_REALLOC(*nbr_sselt_boco_som, cpt_boco_som + 1, ecs_int_t);
      ECS_REALLOC(*tab_sselt_boco_som, cpt_boco_som + 1, ecs_int_t *);

      (*ind_nom_boco_som)[cpt_boco_som] = ptr_boco->ind_nom;
      (*nbr_sselt_boco_som)[cpt_boco_som] = 0;

      cpt_sselt_boco = 0;

      /*
        Boucle inférieure pour parcourir toutes les C.L. de même nom
        (utile notamment si l'on a plusieurs zones)
      */

      for (ind_boco_sub = ind_boco;
           ind_boco_sub < nbr_boco_tot;
           ind_boco_sub++) {

        ptr_boco = tab_boco + ind_boco_sub;

        if (ptr_boco->ind_nom + 1 != num_boco)
          continue;

        if (ptr_boco->support == CS_CG_ENUM(Vertex)) {

          ptr_zone = tab_zone + ptr_boco->num_zone - 1;

          ecs_loc_pre_cgns__marque_som(base_maillage,
                                       num_boco,
                                       ptr_zone,
                                       ptr_boco->ptset_type,
                                       ptr_boco->npnts,
                                       ptr_boco->pnts,
                                       indic_som);

          /* Marquage C.L. après utilisation et libération mémoire partielle */

          if (ind_boco_sub > ind_boco)
            ptr_boco->ind_nom = - 1;

          ptr_boco->ptset_type = CS_CG_ENUM(PointSetTypeNull);
          ptr_boco->npnts = 0;
          ECS_FREE(ptr_boco->pnts);

        }

      } /* Fin de la boucle inférieure */

      /* On marque les éléments inférieurs effectivement utilisés */

      for (ind_sselt = nbr_sselt_ini;
           ind_sselt < nbr_elt_ent[ECS_ENTMAIL_FAC];
           ind_sselt++) {

        bool_cree = true;
        for (ind_pos_sselt = pos_som_sselt[ind_sselt    ] - 1;
             ind_pos_sselt < pos_som_sselt[ind_sselt + 1] - 1;
             ind_pos_sselt++) {
          if (indic_som[val_som_sselt[ind_pos_sselt] - 1] != num_boco)
            bool_cree = false;
        }

        if (bool_cree == true) {

          if (indic_sselt[ind_sselt] == 0) {
            indic_sselt[ind_sselt] = 1;
            cpt_sselt++;
          }

          liste_sselt[cpt_sselt_boco] = ind_sselt;
          cpt_sselt_boco++;

        }

      }

      if (cpt_sselt_boco > 0) {

        (*nbr_sselt_boco_som)[cpt_boco_som] = cpt_sselt_boco;

        ECS_MALLOC((*tab_sselt_boco_som)[cpt_boco_som],
                   cpt_sselt_boco, ecs_int_t);

        for (ind_sselt = 0; ind_sselt < cpt_sselt_boco; ind_sselt++)
          (*tab_sselt_boco_som)[cpt_boco_som][ind_sselt]
            = liste_sselt[ind_sselt];

        cpt_boco_som++;

      }

    }

  }

  *nbr_boco_som = cpt_boco_som;

  ECS_FREE(indic_som);
  ECS_FREE(liste_sselt);


  /* Compactage des définitions des faces */
  /*--------------------------------------*/

  cpt_sselt      = nbr_sselt_ini;
  cpt_val_sselt  = nbr_val_sselt_ini;

  ECS_MALLOC(renum_sselt, nbr_elt_ent[ECS_ENTMAIL_FAC], ecs_int_t);

  for (ind_sselt = nbr_sselt_ini;
       ind_sselt < nbr_elt_ent[ECS_ENTMAIL_FAC];
       ind_sselt++) {

    if (indic_sselt[ind_sselt] != 0) {

      renum_sselt[ind_sselt] = cpt_sselt;

      pos_som_sselt[cpt_sselt] = cpt_val_sselt + 1;

      for (ind_pos_sselt = pos_som_sselt[ind_sselt    ] - 1;
           ind_pos_sselt < pos_som_sselt[ind_sselt + 1] - 1;
           ind_pos_sselt++) {
        val_som_sselt[cpt_val_sselt++]
          = val_som_sselt[ind_pos_sselt];
      }

      cpt_sselt += 1;

    }
    else

      renum_sselt[ind_sselt] = -1;

  }

  pos_som_sselt[cpt_sselt] = cpt_val_sselt + 1;

  /* Renumérotation des listes */

  for (ind_boco = 0; ind_boco < *nbr_boco_som; ind_boco++) {

    for (ind_sselt = 0;
         ind_sselt < (size_t)((*nbr_sselt_boco_som)[ind_boco]);
         ind_sselt++)

      (*tab_sselt_boco_som)[ind_boco][ind_sselt]
        = renum_sselt[(*tab_sselt_boco_som)[ind_boco][ind_sselt]];

  }

  ECS_FREE(indic_sselt);
  ECS_FREE(renum_sselt);

  ECS_REALLOC(elt_pos_som_ent[ECS_ENTMAIL_FAC], cpt_sselt + 1, ecs_size_t);
  ECS_REALLOC(elt_val_som_ent[ECS_ENTMAIL_FAC], cpt_val_sselt, ecs_int_t);

  nbr_elt_ent[ECS_ENTMAIL_FAC] = cpt_sselt;
}

/*----------------------------------------------------------------------------
 *                        Lecture des éléments
 *----------------------------------------------------------------------------*/

static void
ecs_loc_pre_cgns__init_elt_type(ecs_cgns_elt_t elt_type[NofValidElementTypes])
{
  int i, j;

  /* Initialize all types to ECS_ELT_TYP_NUL so as to handle defaults;
     MIXED will be handled separately in caller code. */

  for (i = 0; i < NofValidElementTypes; i++) {
    elt_type[i].ecs_type = ECS_ELT_TYP_NUL;
    elt_type[i].nbr_som = 0;
    for (j = 0; j < 8; j++)
      elt_type[i].num_som[j] = -1;
  }

  /* Map known and supported element types */

  elt_type[CS_CG_ENUM(ElementTypeNull)].nbr_som = 1;
  elt_type[CS_CG_ENUM(ElementTypeUserDefined)].nbr_som = 1;
  elt_type[CS_CG_ENUM(NODE)].nbr_som = 1;
  elt_type[CS_CG_ENUM(BAR_2)].nbr_som = 2;
  elt_type[CS_CG_ENUM(BAR_3)].nbr_som = 3;

  elt_type[CS_CG_ENUM(TRI_3)].ecs_type = ECS_ELT_TYP_FAC_TRIA;
  elt_type[CS_CG_ENUM(TRI_3)].nbr_som = 3;

  elt_type[CS_CG_ENUM(TRI_6)].ecs_type = ECS_ELT_TYP_FAC_TRIA;
  elt_type[CS_CG_ENUM(TRI_6)].nbr_som = 6;

  elt_type[CS_CG_ENUM(QUAD_4)].ecs_type = ECS_ELT_TYP_FAC_QUAD;
  elt_type[CS_CG_ENUM(QUAD_4)].nbr_som = 4;

  elt_type[CS_CG_ENUM(QUAD_8)].ecs_type = ECS_ELT_TYP_FAC_QUAD;
  elt_type[CS_CG_ENUM(QUAD_8)].nbr_som = 8;

  elt_type[CS_CG_ENUM(QUAD_9)].ecs_type = ECS_ELT_TYP_FAC_QUAD;
  elt_type[CS_CG_ENUM(QUAD_9)].nbr_som = 9;

  elt_type[CS_CG_ENUM(TETRA_4)].ecs_type = ECS_ELT_TYP_CEL_TETRA;
  elt_type[CS_CG_ENUM(TETRA_4)].nbr_som = 4;

  elt_type[CS_CG_ENUM(TETRA_10)].ecs_type = ECS_ELT_TYP_CEL_TETRA;
  elt_type[CS_CG_ENUM(TETRA_10)].nbr_som = 10;

  elt_type[CS_CG_ENUM(PYRA_5)].ecs_type = ECS_ELT_TYP_CEL_PYRAM;
  elt_type[CS_CG_ENUM(PYRA_5)].nbr_som = 5;

  elt_type[CS_CG_ENUM(PYRA_14)].ecs_type = ECS_ELT_TYP_CEL_PYRAM;
  elt_type[CS_CG_ENUM(PYRA_14)].nbr_som = 14;

  elt_type[CS_CG_ENUM(PENTA_6)].ecs_type = ECS_ELT_TYP_CEL_PRISM;
  elt_type[CS_CG_ENUM(PENTA_6)].nbr_som = 6;

  elt_type[CS_CG_ENUM(PENTA_15)].ecs_type = ECS_ELT_TYP_CEL_PRISM;
  elt_type[CS_CG_ENUM(PENTA_15)].nbr_som = 15;

  elt_type[CS_CG_ENUM(PENTA_18)].ecs_type = ECS_ELT_TYP_CEL_PRISM;
  elt_type[CS_CG_ENUM(PENTA_18)].nbr_som = 18;

  elt_type[CS_CG_ENUM(HEXA_8)].ecs_type = ECS_ELT_TYP_CEL_HEXA;
  elt_type[CS_CG_ENUM(HEXA_8)].nbr_som = 8;

  elt_type[CS_CG_ENUM(HEXA_20)].ecs_type = ECS_ELT_TYP_CEL_HEXA;
  elt_type[CS_CG_ENUM(HEXA_20)].nbr_som = 20;

  elt_type[CS_CG_ENUM(HEXA_27)].ecs_type = ECS_ELT_TYP_CEL_HEXA;
  elt_type[CS_CG_ENUM(HEXA_27)].nbr_som = 27;

#if (CGNS_VERSION >= 3000)
  elt_type[CS_CG_ENUM(PYRA_13)].ecs_type = ECS_ELT_TYP_CEL_PYRAM;
  elt_type[CS_CG_ENUM(PYRA_13)].nbr_som = 13;
#endif

  elt_type[CS_CG_ENUM(NGON_n)].ecs_type = ECS_ELT_TYP_FAC_POLY;

#if (CGNS_VERSION >= 3000)
  elt_type[CS_CG_ENUM(NFACE_n)].ecs_type = ECS_ELT_TYP_CEL_POLY;
#endif

  /* Elements vertices in current CGNS versions have the same
     local numbering as that of Code_Saturne. */

  for (i = 0; i < NofValidElementTypes; i++) {
    int nbr_som = ecs_fic_elt_typ_liste_c[elt_type[i].ecs_type].nbr_som;
    for (j = 0; j < nbr_som; j++)
      elt_type[i].num_som[j] = j+1;
  }

}

/*----------------------------------------------------------------------------
 *                        Lecture des éléments
 *----------------------------------------------------------------------------*/

static void
ecs_loc_pre_cgns__lit_ele(ecs_maillage_t             *maillage,
                          const ecs_loc_cgns_base_t  *base_maillage,
                          int                         nzones,
                          ecs_loc_cgns_zone_t        *tab_zone,
                          int                         nbr_boco_tot,
                          ecs_loc_cgns_boco_t        *tab_boco,
                          ecs_int_t                  *nbr_boco_som,
                          ecs_int_t                 **ind_nom_boco_som,
                          ecs_int_t                 **nbr_sselt_boco_som,
                          ecs_int_t                ***tab_sselt_boco_som,
                          ecs_int_t                 **ind_section_cel,
                          ecs_int_t                 **ind_section_fac,
                          ecs_int_t                 **ind_zone_cel,
                          ecs_int_t                 **ind_zone_fac)
{
  ecs_int_t    ient;
  ecs_int_t    ient_max;

  ecs_int_t    ind_type = -1;
  ecs_int_t    ind_pos;
  ecs_int_t    ind_val;
  ecs_int_t    ind_som;
  ecs_int_t    cpt_elt_loc;
  ecs_int_t    cpt_section;
  ecs_int_t    nbr_elt_loc;
  ecs_int_t    nbr_elt_zone;
  ecs_int_t    nbr_som_elt = 0;

  ecs_int_t    ii, jj, kk;
  ecs_int_t    ni, nj, nk;

  ecs_int_t    num_som_deb;
  ecs_int_t    num_som_loc;

  ecs_elt_typ_t  ecs_typ;
  ecs_cgns_elt_t ecs_cgns_elt_liste_c[NofValidElementTypes];

  /* Déclarations des variables de stockage        */
  /* avant transfert dans la structure du maillage */
  /*-----------------------------------------------*/

  size_t       cpt_elt_ent[ECS_N_ENTMAIL];
  size_t       cpt_val_ent[ECS_N_ENTMAIL];
  ecs_int_t    cpt_coul_ent[ECS_N_ENTMAIL];

  ecs_size_t  * elt_pos_som_ent    [ECS_N_ENTMAIL];  /* Positions sommets    */
  ecs_int_t   * elt_val_som_ent    [ECS_N_ENTMAIL];  /* Numéros des sommets  */
  ecs_int_t   * elt_val_couleur_ent[ECS_N_ENTMAIL];  /* Couleurs éléments    */
  ecs_int_t   * val_coul_ent       [ECS_N_ENTMAIL];
  ecs_size_t  * cpt_elt_coul_ent   [ECS_N_ENTMAIL];
  ecs_int_t   * ind_zone_ent       [ECS_N_ENTMAIL];
  ecs_int_t   * ind_section_ent    [ECS_N_ENTMAIL];

  ecs_loc_cgns_section_t  *ptr_section;
  ecs_loc_cgns_zone_t     *ptr_zone;

  /* Déclarations des variables pour CGNS */
  /*-------------------------------------*/

  int         cel_dim;
  int         ind_zone;
  int         ind_section;
  int         num_fic;
  int         num_base;
  int         num_section;
  int         num_zone;
  cgsize_t    taille_loc;

  cgsize_t   *ptr_ele;
  cgsize_t   *parentdata;

  int         ret = 0;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Initialisations */

  num_fic  = base_maillage->num_fic;
  num_base = base_maillage->num_base;

  cel_dim = base_maillage->dim_entite;

  if (cel_dim == 1)
    ient_max = ECS_ENTMAIL_NONE;
  else if (cel_dim == 2)
    ient_max = ECS_ENTMAIL_FAC;
  else
    ient_max = ECS_ENTMAIL_CEL;

  ecs_loc_pre_cgns__init_elt_type(ecs_cgns_elt_liste_c);

  /* Allocations des tableaux locaux */

  /* Attention au décalage de `1' !!!         */
  /* On n'allouera pas les tableaux locaux pour */
  /* `ECS_ENTMAIL_DEB = ECS_ENTMAIL_SOM'      */

  for (ient = 0; ient < ECS_N_ENTMAIL; ient++) {

    cpt_elt_ent        [ient] = 0;
    cpt_val_ent        [ient] = 0;
    cpt_coul_ent       [ient] = 0;

    elt_pos_som_ent    [ient] = NULL;
    elt_val_som_ent    [ient] = NULL;
    elt_val_couleur_ent[ient] = NULL;

    val_coul_ent       [ient] = NULL;
    cpt_elt_coul_ent   [ient] = NULL;

    ind_zone_ent       [ient] = NULL;
    ind_section_ent    [ient] = NULL;

  }

  /*---------------------------*/
  /* Lecture des connectivités */
  /*---------------------------*/

  for (ind_zone = 0; ind_zone < nzones; ind_zone++) {

    num_zone = ind_zone + 1;
    ptr_zone = tab_zone + ind_zone;

    for (ind_section = 0;
         ind_section < ptr_zone->nbr_sections;
         ind_section++) {

      num_section = ind_section + 1;
      ptr_section = ptr_zone->tab_sections + ind_section;

      ret = cg_ElementDataSize(num_fic, num_base, num_zone, num_section,
                               &taille_loc);

      /*
        Si la section est de type "MIXED", taille_loc correspond
        aux tailles "pos" + "val"; sinon, taille_loc correspond
        aux tailles "val" uniquement. On ne sait pas si ces
        données correspondent aux éléments de dimension maximale,
        et on les lit donc immédiatement.
      */


      if (ret == CG_OK) {

        ECS_MALLOC(ptr_section->elems, taille_loc, cgsize_t);

        if (ptr_section->parent > 0)
          ECS_MALLOC(parentdata,
                     (  ptr_section->num_elt_fin
                      - ptr_section->num_elt_deb + 1) * 4,
                     cgsize_t);
        else
          parentdata = NULL;

#if CGNS_VERSION >= 3400

        if (   ptr_section->type == CS_CG_ENUM(MIXED)
            || ptr_section->type == CS_CG_ENUM(NFACE_n)
            || ptr_section->type == CS_CG_ENUM(NGON_n)) {
          nbr_elt_loc = ptr_section->num_elt_fin - ptr_section->num_elt_deb + 1;
          ECS_MALLOC(ptr_section->offsets, nbr_elt_loc + 1, cgsize_t);
          ret = cg_poly_elements_read(num_fic, num_base, num_zone, num_section,
                                      ptr_section->elems, ptr_section->offsets,
                                      parentdata);
        }
        else
          ret = cg_elements_read(num_fic, num_base, num_zone, num_section,
                                 ptr_section->elems, parentdata);

#else

        ret = cg_elements_read(num_fic, num_base, num_zone, num_section,
                               ptr_section->elems, parentdata);

#endif

        if (ret != CG_OK)
          ecs_error(__FILE__, __LINE__, 0,
                    _("CGNS error:\n%s"), cg_get_error());

        ECS_FREE(parentdata);

      }

    }

  }

  /*--------------------------------------------------------*/
  /* Dimensionnement de la connectivité nodale des éléments */
  /*--------------------------------------------------------*/

  for (ind_zone = 0; ind_zone < nzones; ind_zone++) {

    num_zone = ind_zone + 1;
    ptr_zone = tab_zone + ind_zone;

    nbr_elt_zone = 0;

    cgsize_t max_elt_zone = 0;

    /* Cas d'une zone structurée */

    if (ptr_zone->type == CS_CG_ENUM(Structured)) {

      nbr_elt_loc   = ptr_zone->nbr_cel;
      nbr_elt_zone += nbr_elt_loc;

      if (cel_dim == 2) {
        cpt_elt_ent[ECS_ENTMAIL_FAC] += nbr_elt_loc;
        cpt_val_ent[ECS_ENTMAIL_FAC] += nbr_elt_loc * 4;
      }
      else if (cel_dim == 3) {
        cpt_elt_ent[ECS_ENTMAIL_CEL] += nbr_elt_loc;
        cpt_val_ent[ECS_ENTMAIL_CEL] += nbr_elt_loc * 8;
      }

      max_elt_zone = nbr_elt_zone;

    }

    /* Cas d'une zone non structurée */

    else if (ptr_zone->type == CS_CG_ENUM(Unstructured)) {

      for (ind_section = 0;
           ind_section < ptr_zone->nbr_sections;
           ind_section++) {

        num_section = ind_section + 1;
        ptr_section = ptr_zone->tab_sections + ind_section;

        nbr_elt_loc   = ptr_section->num_elt_fin - ptr_section->num_elt_deb + 1;
        nbr_elt_zone += nbr_elt_loc;

        ecs_typ     = ecs_cgns_elt_liste_c[ptr_section->type].ecs_type;
        nbr_som_elt = ecs_fic_elt_typ_liste_c[ecs_typ].nbr_som;

        ient = ecs_maillage_pre__ret_typ_geo(ecs_typ);

        max_elt_zone = ECS_MAX(max_elt_zone, ptr_section->num_elt_fin);

        if (ptr_section->type == CS_CG_ENUM(MIXED)) {

          cpt_elt_loc = 0;

          cgsize_t *elt_idx = ptr_section->offsets;
          ptr_ele = ptr_section->elems;

          while (cpt_elt_loc < nbr_elt_loc) {

            ind_type    = *ptr_ele;

#if (CGNS_VERSION < 3200)
            if (ind_type < CS_CG_ENUM(NGON_n)) {
              ecs_typ     = ecs_cgns_elt_liste_c[ind_type].ecs_type;
              nbr_som_elt = ecs_fic_elt_typ_liste_c[ecs_typ].nbr_som;
              ptr_ele += ecs_cgns_elt_liste_c[ind_type].nbr_som + 1;
            }
            else {
              ecs_typ     = ECS_ELT_TYP_FAC_POLY;
              nbr_som_elt = *ptr_ele - CS_CG_ENUM(NGON_n);
              ptr_ele += nbr_som_elt + 1;
            }
#else
            ecs_typ     = ecs_cgns_elt_liste_c[ind_type].ecs_type;
            nbr_som_elt = ecs_fic_elt_typ_liste_c[ecs_typ].nbr_som;
            if (ind_type != CS_CG_ENUM(NGON_n))
              ptr_ele += ecs_cgns_elt_liste_c[ind_type].nbr_som + 1;
            else if (elt_idx != NULL) {
              nbr_som_elt = elt_idx[cpt_elt_loc+1] - elt_idx[cpt_elt_loc+1] - 1;
              ptr_ele += nbr_som_elt + 1;
            }
            else
              ecs_error(__FILE__, __LINE__, 0,
                        _("CGNS: unhandled NGON_n element in MIXED section\n"));
#endif

            cpt_elt_loc++;

            /*
              Calcul entité correspondante sans passer par
              ecs_entmail_pre__ret_typ_geo(ecs_typ) pour éviter
              appel de fonction à chaque élément.
            */

            if (ecs_typ == ECS_ELT_TYP_NUL)
              ient = ECS_ENTMAIL_NONE;
            else if (ecs_typ < ECS_ELT_TYP_CEL_TETRA)
              ient = ECS_ENTMAIL_FAC;
            else if (ecs_typ < ECS_ELT_TYP_FAC_POLY)
              ient = ECS_ENTMAIL_CEL;
            else if (ecs_typ == ECS_ELT_TYP_FAC_POLY)
              ient = ECS_ENTMAIL_FAC;
            else if (ecs_typ == ECS_ELT_TYP_CEL_POLY)
              ient = ECS_ENTMAIL_CEL;
            else
              ient = ECS_ENTMAIL_NONE;

            if (ient != ECS_ENTMAIL_NONE) {
              cpt_elt_ent[ient] += 1;
              cpt_val_ent[ient] += nbr_som_elt;
            }

          }
        }

        else if (ptr_section->type == CS_CG_ENUM(NGON_n)) {

#if CGNS_VERSION >= 3400

          cpt_elt_ent[ient] += nbr_elt_loc;
          cpt_val_ent[ient] += ptr_section->offsets[nbr_elt_loc];

#else

          cpt_elt_loc = 0;
          ptr_ele = ptr_section->elems;

          while (cpt_elt_loc < nbr_elt_loc) {

            nbr_som_elt = *ptr_ele;
            ptr_ele += nbr_som_elt + 1;

            cpt_elt_loc++;

            cpt_elt_ent[ient] += 1;
            cpt_val_ent[ient] += nbr_som_elt;

          }

#endif

        }

#if CGNS_VERSION >= 3000

        else if (ptr_section->type == CS_CG_ENUM(NFACE_n)) {

#if CGNS_VERSION >= 3400

          cpt_elt_ent[ient] += nbr_elt_loc;

#else

          cpt_elt_loc = 0;
          ptr_ele = ptr_section->elems;

          while (cpt_elt_loc < nbr_elt_loc) {

            ecs_int_t nbr_fac_elt = *ptr_ele;
            ptr_ele += nbr_fac_elt + 1;

            cpt_elt_loc++;

            cpt_elt_ent[ient] += 1;

            /* cpt_val_ent[ient] handled later */

          }

#endif

        }

#endif

        else { /* Regular section */

          if (ient != ECS_ENTMAIL_NONE) {
            cpt_elt_ent[ient] += nbr_elt_loc;
            cpt_val_ent[ient] += nbr_elt_loc * nbr_som_elt;
          }

        }

      } /* Fin boucle sur les sections */

    } /* Fin traitement structuré/non structuré */

    /* Ajout pour conditions aux limites */

    if (ptr_zone->renum_size > -1) {
      cgsize_t n = max_elt_zone + 1;
      ptr_zone->renum_size = n;
      ECS_MALLOC(ptr_zone->renum, n, cgsize_t);
      for (cgsize_t i = 0; i < n; i++)
        ptr_zone->renum[i] = 0;
    }

  } /* Fin boucle sur les zones */


  /* Allocation mémoire et remise compteurs à zéro */
  /*-----------------------------------------------*/

  for (ient = 0; ient < ECS_N_ENTMAIL; ient++) {

    if (cpt_elt_ent[ient] > 0) {

      ECS_MALLOC(elt_pos_som_ent[ient], cpt_elt_ent[ient] + 1, ecs_size_t);
      ECS_MALLOC(elt_val_som_ent[ient], cpt_val_ent[ient], ecs_int_t);

      if (ient == ient_max - 1) {
        if (ind_zone_fac != NULL) {
          ECS_MALLOC(ind_zone_ent[ient], cpt_elt_ent[ient], ecs_int_t);
          *ind_zone_fac = ind_zone_ent[ient];
        }
        if (ind_section_fac != NULL) {
          ECS_MALLOC(ind_section_ent[ient], cpt_elt_ent[ient], ecs_int_t);
          *ind_section_fac = ind_section_ent[ient];
        }
      }
      else if (ient == ient_max) {
        if (ind_zone_cel != NULL) {
          ECS_MALLOC(ind_zone_ent[ient], cpt_elt_ent[ient], ecs_int_t);
          *ind_zone_cel = ind_zone_ent[ient];
        }
        if (ind_section_cel != NULL) {
          ECS_MALLOC(ind_section_ent[ient], cpt_elt_ent[ient], ecs_int_t);
          *ind_section_cel = ind_section_ent[ient];
        }
      }

      cpt_elt_ent [ient] = 0;
      cpt_val_ent [ient] = 0;

      elt_pos_som_ent[ient][0] = 1;

    }

  }


  /*-----------------------------------------------------*/
  /* Construction de la connectivité nodale des éléments */
  /*-----------------------------------------------------*/

  cpt_section = 0;
  num_som_deb = 1;

  ecs_int_t face_id_shift = 0;

  for (ind_zone = 0; ind_zone < nzones; ind_zone++) {

    num_zone = ind_zone + 1;
    ptr_zone = tab_zone + ind_zone;

    /* Cas d'une zone structurée */

    if (ptr_zone->type == CS_CG_ENUM(Structured)) {

      nbr_elt_loc = ptr_zone->nbr_cel;
      cpt_elt_loc = 0;

      if (cel_dim == 2) {

        ient = ECS_ENTMAIL_FAC;
        ecs_typ = ECS_ELT_TYP_FAC_QUAD;

        nbr_som_elt = 4;

        ni = ptr_zone->taille[0];
        nj = ptr_zone->taille[1];

        for (jj = 0; jj < nj - 1; jj++) {
          for (ii = 0; ii < ni - 1; ii++) {

            ind_pos = cpt_elt_ent[ient];
            ind_val = elt_pos_som_ent[ient][ind_pos] - 1;
            elt_pos_som_ent[ient][ind_pos + 1]
              =  elt_pos_som_ent[ient][ind_pos] + nbr_som_elt;

            num_som_loc = num_som_deb + ii + jj*ni;

            elt_val_som_ent[ient][ind_val++] = num_som_loc;
            elt_val_som_ent[ient][ind_val++] = num_som_loc + 1;
            elt_val_som_ent[ient][ind_val++] = num_som_loc + 1 + ni;
            elt_val_som_ent[ient][ind_val++] = num_som_loc + ni;

            if (ind_zone_ent[ient] != NULL)
              ind_zone_ent[ient][cpt_elt_ent[ient]] = ind_zone;

            cpt_elt_ent[ient]++;

          }
        }

      }
      else if (cel_dim == 3) {

        ient = ECS_ENTMAIL_CEL;
        ecs_typ = ECS_ELT_TYP_CEL_HEXA;
        nbr_som_elt = 8;

        ni = ptr_zone->taille[0];
        nj = ptr_zone->taille[1];
        nk = ptr_zone->taille[2];

        for (kk = 0; kk < nk - 1; kk++) {
          for (jj = 0; jj < nj - 1; jj++) {
            for (ii = 0; ii < ni - 1; ii++) {

              ind_pos = cpt_elt_ent[ient];
              ind_val = elt_pos_som_ent[ient][ind_pos] - 1;
              elt_pos_som_ent[ient][ind_pos + 1]
                =  elt_pos_som_ent[ient][ind_pos] + nbr_som_elt;

              num_som_loc = num_som_deb + ii + jj*ni + kk*ni*nj;

              elt_val_som_ent[ient][ind_val++] = num_som_loc;
              elt_val_som_ent[ient][ind_val++] = num_som_loc + 1;
              elt_val_som_ent[ient][ind_val++] = num_som_loc + 1 + ni;
              elt_val_som_ent[ient][ind_val++] = num_som_loc + ni;
              elt_val_som_ent[ient][ind_val++] = num_som_loc + ni*nj;
              elt_val_som_ent[ient][ind_val++] = num_som_loc + ni*nj + 1;
              elt_val_som_ent[ient][ind_val++] = num_som_loc + ni*nj + 1 + ni;
              elt_val_som_ent[ient][ind_val++] = num_som_loc + ni*nj + ni;

              if (ind_zone_ent[ient] != NULL)
                ind_zone_ent[ient][cpt_elt_ent[ient]] = ind_zone;
              if (ind_section_ent[ient] != NULL)
                ind_section_ent[ient][cpt_elt_ent[ient]] = 0;

              cpt_elt_ent[ient]++;

            }
          }
        }

      }

    }

    /* Cas d'une zone non structurée */

    else if (ptr_zone->type == CS_CG_ENUM(Unstructured)) {

      for (ind_section = 0;
           ind_section < ptr_zone->nbr_sections;
           ind_section++) {

        num_section = ind_section + 1;
        ptr_section = ptr_zone->tab_sections + ind_section;

        nbr_elt_loc = ptr_section->num_elt_fin - ptr_section->num_elt_deb + 1;
        cpt_elt_loc = 0;

        /* Counting stage for polyhedra */

#if CGNS_VERSION >= 3000

        if (ptr_section->type == CS_CG_ENUM(NFACE_n)) {

          ient = ECS_ENTMAIL_CEL;

#if CGNS_VERSION >= 3400
          cgsize_t *elt_idx = ptr_section->offsets;
#endif
          ptr_ele = ptr_section->elems;

          size_t   connect_size = 0;
          const size_t  *pos_som_fac = elt_pos_som_ent[ECS_ENTMAIL_FAC];

          ecs_int_t n_elts_loc =    ptr_section->num_elt_fin
                                  - ptr_section->num_elt_deb + 1;

          for (ecs_int_t ielt = 0; ielt < n_elts_loc; ielt++) {

#if CGNS_VERSION >= 3400
            ecs_int_t nbr_fac_elt = elt_idx[ielt+1] - elt_idx[ielt];
#else
            ecs_int_t nbr_fac_elt = *ptr_ele;
            ptr_ele += 1;
#endif

            for (ecs_int_t ind_fac = 0; ind_fac < nbr_fac_elt; ind_fac++) {
              ecs_int_t num_fac = *(ptr_ele + ind_fac);
              if (num_fac < 0)
                num_fac = -num_fac;
              num_fac += face_id_shift;
              connect_size +=   pos_som_fac[num_fac]
                              - pos_som_fac[num_fac - 1] + 1;
            }

            ptr_ele += nbr_fac_elt;

          }

          cpt_val_ent[ient] += connect_size;
          ECS_REALLOC(elt_val_som_ent[ient], cpt_val_ent[ient], ecs_int_t);

        }

#endif

        cgsize_t *elt_idx = ptr_section->offsets;

        ptr_ele = ptr_section->elems;

        /* General data for section;
           will be updated for MIXED, NGON_n, and NFACE_n */

        ind_type    = ptr_section->type;
        ecs_typ     = ecs_cgns_elt_liste_c[ind_type].ecs_type;
        nbr_som_elt = ecs_fic_elt_typ_liste_c[ecs_typ].nbr_som;

        ient = ecs_maillage_pre__ret_typ_geo(ecs_typ);

        while (cpt_elt_loc < nbr_elt_loc) {

          if (ptr_section->type == CS_CG_ENUM(MIXED)) {

            ind_type    = *ptr_ele;

#if (CGNS_VERSION < 3200)
            if (ind_type < CS_CG_ENUM(NGON_n)) {
              ecs_typ     = ecs_cgns_elt_liste_c[ind_type].ecs_type;
              nbr_som_elt = ecs_fic_elt_typ_liste_c[ecs_typ].nbr_som;
            }
            else {
              ecs_typ     = ECS_ELT_TYP_FAC_POLY;
              nbr_som_elt = *ptr_ele - CS_CG_ENUM(NGON_n);
            }
#else
            ecs_typ     = ecs_cgns_elt_liste_c[ind_type].ecs_type;
            if (ind_type != CS_CG_ENUM(NGON_n))
              nbr_som_elt = ecs_fic_elt_typ_liste_c[ecs_typ].nbr_som;
            else
              nbr_som_elt = elt_idx[cpt_elt_loc+1] - elt_idx[cpt_elt_loc+1] - 1;
#endif

            ptr_ele += 1;

            /*
              Calcul entité correspondante sans passer par
              ecs_entmail_pre__ret_typ_geo(ecs_typ) pour éviter
              appel de fonction à chaque élément.
            */

            if (ecs_typ == ECS_ELT_TYP_NUL)
              ient = ECS_ENTMAIL_NONE;
            else if (ecs_typ < ECS_ELT_TYP_CEL_TETRA)
              ient = ECS_ENTMAIL_FAC;
            else if (ecs_typ < ECS_ELT_TYP_FAC_POLY)
              ient = ECS_ENTMAIL_CEL;
            else if (ecs_typ == ECS_ELT_TYP_FAC_POLY)
              ient = ECS_ENTMAIL_FAC;
            else if (ecs_typ == ECS_ELT_TYP_CEL_POLY)
              ient = ECS_ENTMAIL_CEL;
            else
              ient = ECS_ENTMAIL_NONE;

            if (ient != ECS_ENTMAIL_NONE) {

              ind_pos = cpt_elt_ent[ient];
              ind_val = elt_pos_som_ent[ient][ind_pos] - 1;
              elt_pos_som_ent[ient][ind_pos + 1]
                =  elt_pos_som_ent[ient][ind_pos] + nbr_som_elt;

#if (CGNS_VERSION < 3200)
              if (ind_type < CS_CG_ENUM(NGON_n)) {
                for (ind_som = 0; ind_som < nbr_som_elt; ind_som++)
                  elt_val_som_ent[ient][ind_val++]
                    = *(ptr_ele
                        + ecs_cgns_elt_liste_c[ind_type].num_som[ind_som] - 1)
                        + num_som_deb - 1;
                ptr_ele += ecs_cgns_elt_liste_c[ind_type].nbr_som;
              }
              else {
                for (ind_som = 0; ind_som < nbr_som_elt; ind_som++)
                  elt_val_som_ent[ient][ind_val++]
                    = *(ptr_ele + ind_som) + num_som_deb - 1;
                ptr_ele += nbr_som_elt;
              }
#else
              for (ind_som = 0; ind_som < nbr_som_elt; ind_som++)
                elt_val_som_ent[ient][ind_val++]
                  = *(ptr_ele
                      + ecs_cgns_elt_liste_c[ind_type].num_som[ind_som] - 1)
                      + num_som_deb - 1;
              ptr_ele += ecs_cgns_elt_liste_c[ind_type].nbr_som;
#endif
            }

          }

          else if (ptr_section->type == CS_CG_ENUM(NGON_n)) {

#if CGNS_VERSION >= 3400
            nbr_som_elt = elt_idx[cpt_elt_loc + 1] - elt_idx[cpt_elt_loc];
#else
            nbr_som_elt = *ptr_ele;
            ptr_ele += 1;
#endif

            ind_pos = cpt_elt_ent[ient];
            ind_val = elt_pos_som_ent[ient][ind_pos] - 1;
            elt_pos_som_ent[ient][ind_pos + 1]
              =  elt_pos_som_ent[ient][ind_pos] + nbr_som_elt;

            for (ind_som = 0; ind_som < nbr_som_elt; ind_som++) {
              elt_val_som_ent[ient][ind_val++]
                = *(ptr_ele + ind_som) + num_som_deb - 1;
            }

            ptr_ele += nbr_som_elt;

          }

#if CGNS_VERSION >= 3000

          else if (ptr_section->type == CS_CG_ENUM(NFACE_n)) {

#if CGNS_VERSION >= 3400
            ecs_int_t nbr_fac_elt =    elt_idx[cpt_elt_loc + 1]
                                     - elt_idx[cpt_elt_loc];
#else
            ecs_int_t nbr_fac_elt = *ptr_ele;
            ptr_ele += 1;
#endif

            ind_pos = cpt_elt_ent[ient];
            ind_val = elt_pos_som_ent[ient][ind_pos] - 1;

            const size_t     *pos_som_fac = elt_pos_som_ent[ECS_ENTMAIL_FAC];
            const ecs_int_t  *val_som_fac = elt_val_som_ent[ECS_ENTMAIL_FAC];

            for (ecs_int_t i = 0; i < nbr_fac_elt; i++) {

              ecs_int_t num_fac = *(ptr_ele + i);
              ecs_int_t ind_fac = ECS_ABS(num_fac) - 1 + face_id_shift;

              size_t s_id = pos_som_fac[ind_fac] - 1;
              size_t e_id = pos_som_fac[ind_fac + 1] - 1;

              if (num_fac < 0) {
                ecs_int_t num_som_deb_fac = val_som_fac[e_id - 1];
                for (size_t j = e_id; j > s_id; j--)
                  elt_val_som_ent[ient][ind_val++] = val_som_fac[j-1];
                elt_val_som_ent[ient][ind_val++] = num_som_deb_fac;
              }

              else {
                ecs_int_t num_som_deb_fac = val_som_fac[s_id];
                for (size_t j = s_id; j < e_id; j++)
                  elt_val_som_ent[ient][ind_val++] = val_som_fac[j];
                elt_val_som_ent[ient][ind_val++] = num_som_deb_fac;
              }

            }

            elt_pos_som_ent[ient][ind_pos + 1] =  ind_val + 1;
            ptr_ele += nbr_fac_elt;

          }

#endif

          else if (ient != ECS_ENTMAIL_NONE) {

            ind_pos = cpt_elt_ent[ient];
            ind_val = elt_pos_som_ent[ient][ind_pos] - 1;
            elt_pos_som_ent[ient][ind_pos + 1]
              =  elt_pos_som_ent[ient][ind_pos] + nbr_som_elt;

            for (ind_som = 0; ind_som < nbr_som_elt; ind_som++) {
              elt_val_som_ent[ient][ind_val++]
                = *(ptr_ele
                    + ecs_cgns_elt_liste_c[ind_type].num_som[ind_som] - 1)
                    + num_som_deb - 1;
            }

            ptr_ele += ecs_cgns_elt_liste_c[ind_type].nbr_som;

          }

          if (ient != ECS_ENTMAIL_NONE) {

            cgsize_t cg_elt_id = ptr_section->num_elt_deb + cpt_elt_loc - 1;

            if (ptr_zone->renum_size > 0) {
              if (ient == ient_max - 1)
                ptr_zone->renum[cg_elt_id] = cpt_elt_ent[ient] + 1;
              else if (ient == ient_max)
                ptr_zone->renum[cg_elt_id] = -(cpt_elt_ent[ient] + 1);
            }

            if (ind_zone_ent[ient] != NULL)
              ind_zone_ent[ient][cpt_elt_ent[ient]] = ind_zone;

            if (ind_section_ent[ient] != NULL)
              ind_section_ent[ient][cpt_elt_ent[ient]] = cpt_section;

            cpt_elt_ent[ient]++;

          }

          cpt_elt_loc++;

        }

        /* On libère la mémoire */

        ECS_FREE(ptr_section->elems);

        cpt_section += 1;

      } /* Fin boucle sur les sections */

    } /* Fin traitement structuré/non structuré */

    num_som_deb += ptr_zone->nbr_som;

    face_id_shift = cpt_elt_ent[ECS_ENTMAIL_FAC];

  } /* Fin boucle sur les zones */


  /* Suppression des entités inutilisées */

  for (ient = 0; ient < ECS_N_ENTMAIL; ient++) {

    if (cpt_elt_ent[ient] > 0 && ient < ient_max - 1) {

      ECS_FREE(elt_pos_som_ent[ient]);
      ECS_FREE(elt_val_som_ent[ient]);

      cpt_elt_ent [ient] = 0;
      cpt_val_ent [ient] = 0;

    }

  }

  /* Création des entités supplémentaires pour porter les C.L. aux sommets */
  /*-----------------------------------------------------------------------*/

  ecs_loc_pre__cgns__cree_ent_inf_som(base_maillage,
                                      cpt_elt_ent,
                                      elt_pos_som_ent,
                                      elt_val_som_ent,
                                      nzones,
                                      tab_zone,
                                      nbr_boco_tot,
                                      tab_boco,
                                      nbr_boco_som,
                                      ind_nom_boco_som,
                                      nbr_sselt_boco_som,
                                      tab_sselt_boco_som);

  /* Transfert des valeurs lues dans les structures d'entité de maillage */
  /*=====================================================================*/

  ecs_maillage_pre__cree_elt(maillage,
                             cpt_elt_ent,
                             elt_pos_som_ent,
                             elt_val_som_ent,
                             NULL,
                             elt_val_couleur_ent,
                             cpt_coul_ent,
                             val_coul_ent,
                             cpt_elt_coul_ent);
}

/*----------------------------------------------------------------------------
 *                   Création des groupes basés sur les C.L.
 *----------------------------------------------------------------------------*/

static void
ecs_loc_pre_cgns__cree_grps_boco(const ecs_loc_cgns_base_t  *base_maillage,
                                 ecs_loc_cgns_zone_t        *tab_zone,
                                 int                         nbr_boco_tot,
                                 char                      **tab_nom_boco,
                                 ecs_loc_cgns_boco_t        *tab_boco,
                                 ecs_int_t                   nbr_boco_som,
                                 ecs_int_t                 **ind_nom_boco_som,
                                 ecs_int_t                 **nbr_sselt_boco_som,
                                 ecs_int_t                ***tab_sselt_boco_som,
                                 size_t                     *ent_nbr_elt,
                                 ecs_table_t               **ent_table_grp)
{
  ecs_int_t    ient;
  ecs_int_t    ient_max;

  ecs_int_t    ind;
  ecs_int_t    ind_boco;
  ecs_int_t    ind_boco_sub;
  ecs_int_t    ind_ent;
  ecs_int_t    ind_glob;
  ecs_int_t    ind_nom;

  int          cel_dim;

  ecs_int_t    nbr_sselt_boco;
  ecs_int_t   *tab_sselt_boco;

  ecs_loc_cgns_zone_t   *ptr_zone;
  ecs_loc_cgns_boco_t   *ptr_boco;

  /* Déclarations des variables de stockage        */
  /* avant transfert dans la structure du maillage */
  /*-----------------------------------------------*/

  char       *nom_boco;

  bool         bool_aff_grp;

  ecs_descr_t  *descr_grp;
  size_t        ent_cpt_elt[ECS_N_ENTMAIL];
  ecs_int_t    *ent_val_grp[ECS_N_ENTMAIL];

  ecs_table_t  *table_grp;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Initialisations */

  cel_dim = base_maillage->dim_entite;

  if (cel_dim == 1)
    return;

  else if (cel_dim == 2)
    ient_max = ECS_ENTMAIL_FAC;
  else
    ient_max = ECS_ENTMAIL_CEL;

  ind_boco = 0;

  /* Allocations des tableaux locaux */

  for (ient = 0; ient < ECS_N_ENTMAIL; ient++) {

    ent_table_grp[ient] = NULL;

    ent_cpt_elt[ient] = 0;
    ent_val_grp[ient] = NULL;

  }

  /*-----------------------------------------*/
  /* Génération des groupes associés aux C.L */
  /*-----------------------------------------*/

  for (ind_boco = 0; ind_boco < nbr_boco_tot; ind_boco++) {

    ptr_boco = tab_boco + ind_boco;

    ind_nom  = ptr_boco->ind_nom;

    if (ind_nom < 0)
      continue;

    nom_boco = tab_nom_boco[ind_nom];

    /* On alloue et initialise pour le groupe à lire */

    for (ient = ient_max - 1; ient <= ient_max; ient++) {

      ent_cpt_elt[ient] = 0;

      if (ent_nbr_elt[ient] != 0)
        ECS_MALLOC(ent_val_grp[ient], ent_nbr_elt[ient], ecs_int_t);

      for (ind = 0; ind < (ecs_int_t)(ent_nbr_elt[ient]); ind++)
        ent_val_grp[ient][ind] = 0;

    }

    /* Boucle pour le traitement des C.L. aux sommets */
    /*------------------------------------------------*/

    for (ind_boco_sub = 0;
         ind_boco_sub < nbr_boco_som;
         ind_boco_sub++) {

      if ((*ind_nom_boco_som)[ind_boco_sub] != ind_nom)
        continue;

      /* Traitement préparé dans ecs_loc_pre__cgns__cree_ent_inf_som() */

      ient = ient_max - 1;

      nbr_sselt_boco = (*nbr_sselt_boco_som)[ind_boco_sub];
      tab_sselt_boco = (*tab_sselt_boco_som)[ind_boco_sub];

      for (ind = 0; ind < nbr_sselt_boco; ind++) {

        ind_glob = tab_sselt_boco[ind];

        /* Stockage des valeurs lues avant transfert dans maillage */

        ent_val_grp[ient][ind_glob] = 1;
        ent_cpt_elt[ient]++;

      }

      nbr_sselt_boco = 0;
      tab_sselt_boco = NULL;

      (*nbr_sselt_boco_som)[ind_boco_sub] = 0;
      ECS_FREE((*tab_sselt_boco_som)[ind_boco_sub]);

      /* Marquage C.L. sommets comme déjà utilisée */

      (*ind_nom_boco_som)[ind_boco_sub] = -1;

    }

    /* Boucle pour le traitement des C.L. aux faces */
    /*----------------------------------------------*/

    /*
      Boucle inférieure pour parcourir toutes les C.L. de même nom
      (utile notamment si l'on a plusieurs zones)
    */

    for (ind_boco_sub = ind_boco;
         ind_boco_sub < nbr_boco_tot;
         ind_boco_sub++) {

      ptr_boco = tab_boco + ind_boco_sub;

      if (ptr_boco->ind_nom != ind_nom)
        continue;

      /* Traitement selon le support */
      /*-----------------------------*/

      if (ptr_boco->support == CS_CG_ENUM(FaceCenter)) {

        ptr_zone = tab_zone + ptr_boco->num_zone - 1;

        /* Liste définie par numéro de début et fin */

        if (   ptr_boco->ptset_type == CS_CG_ENUM(PointRange)
            || ptr_boco->ptset_type == CS_CG_ENUM(ElementRange)) {

          for (ind = ptr_boco->pnts[0]; ind <= ptr_boco->pnts[1]; ind++) {

            ind_ent = ind - ptr_zone->num_elt_deb;

            /* Stockage des valeurs lues avant transfert dans maillage */

            if (ptr_zone->renum[ind_ent] > 0) {
              ient = ient_max - 1;
              ind_glob = ptr_zone->renum[ind_ent] - 1;
              ent_val_grp[ient][ind_glob] = 1;
              ent_cpt_elt[ient]++;
            }
            else if (ptr_zone->renum[ind_ent] < 0) {
              ient = ient_max;
              ind_glob = - ptr_zone->renum[ind_ent] - 1;
              ent_val_grp[ient][ind_glob] = 1;
              ent_cpt_elt[ient]++;
            }

          }

        }

        /* Liste définie explicitement */

        else {

          for (ind = 0; ind < ptr_boco->npnts; ind++) {

            ind_ent = ptr_boco->pnts[ind] - ptr_zone->num_elt_deb;

            /* If boundary condition references elements not present,
               ignore it (workaround for bug in ICEM Meshing 13 output). */

            if (ind_ent >= ptr_zone->renum_size)
              continue;

            /* Stockage des valeurs lues avant transfert dans maillage */

            if (ptr_zone->renum[ind_ent] > 0) {
              ient = ient_max - 1;
              ind_glob = ptr_zone->renum[ind_ent] - 1;
              ent_val_grp[ient][ind_glob] = 1;
              ent_cpt_elt[ient]++;
            }
            else if (ptr_zone->renum[ind_ent] < 0) {
              ient = ient_max;
              ind_glob = - ptr_zone->renum[ind_ent] - 1;
              ent_val_grp[ient][ind_glob] = 1;
              ent_cpt_elt[ient]++;
            }

          }

        }

        /* Marquage C.L. après utilisation et libération mémoire partielle */

        ptr_boco->ind_nom = - 1;

        ECS_FREE(ptr_boco->pnts);

      }

    } /* Fin de la boucle inférieure */

    /* Retour si aucune C.L. traitée */

    if (ient > ient_max)
      continue;

    /* Remplissage des entités du maillage */
    /*-------------------------------------*/

    bool_aff_grp = false;

    if (ent_cpt_elt[ient] != 0) {

      bool_aff_grp = true;

      assert(ent_cpt_elt[ient] <= (size_t)(ent_nbr_elt[ient]));

      /* Création du descripteur de table correspondant au groupe lu */

      descr_grp = ecs_descr__cree(ECS_DESCR_IDE_NUL,
                                  nom_boco);

      /* Transformation du tableau référencant le groupe en une table */

      table_grp = ecs_table__transforme_tableau(ent_nbr_elt[ient],
                                                ent_val_grp[ient],
                                                descr_grp);

      if (ent_table_grp[ient] != NULL)
          ecs_table_att__assemble(ent_table_grp[ient],
                                  table_grp);

      else
        ent_table_grp[ient] = table_grp;

    } /* Fin si le nombre d'éléments référencant le groupe n'est pas nul */

    /* Libération mémoire */

    for (ient = ient_max - 1; ient <= ient_max; ient++) {
      if (ent_nbr_elt[ient] != 0)
        ECS_FREE(ent_val_grp[ient]);
    }

    /* Affichage du bilan des données lues pour les groupes */
    /*------------------------------------------------------*/

    if (bool_aff_grp == true)
      printf("  %s \"%s\"\n",
             _("Group"), nom_boco);

    ecs_maillage_pre__aff_nbr_par_ent(0,
                                      ent_cpt_elt,
                                      0);

  } /* Fin de la boucle sur les C.L. */

  /* Libération mémoire associée aux tableaux sur les sommets */

  if (nbr_boco_som > 0) {
    ECS_FREE(*ind_nom_boco_som);
    ECS_FREE(*nbr_sselt_boco_som);
    ECS_FREE(*tab_sselt_boco_som);
  }

  /* Libération de la mémoire associée à d'éventuelles listes non utilisées */

  for (ind_boco = 0; ind_boco < nbr_boco_tot; ind_boco++)
    ECS_FREE(tab_boco[ind_boco].pnts);
}

/*----------------------------------------------------------------------------
 * Création des groupes basés sur les zones et sections
 *----------------------------------------------------------------------------*/

static void
ecs_loc_pre_cgns__cree_grps_zs(const ecs_loc_cgns_base_t  *base_maillage,
                               int                         nzones,
                               ecs_loc_cgns_zone_t        *tab_zone,
                               ecs_int_t                  *tab_ind_section_cel,
                               ecs_int_t                  *tab_ind_section_fac,
                               ecs_int_t                  *tab_ind_zone_cel,
                               ecs_int_t                  *tab_ind_zone_fac,
                               size_t                     *ent_nbr_elt,
                               ecs_table_t               **ent_table_grp)
{
  int          ient;
  int          ient_max;

  ecs_int_t    ind;

  int          cel_dim;

  ecs_loc_cgns_zone_t      *ptr_zone;
  ecs_loc_cgns_section_t   *ptr_section;

  ecs_int_t     nbr_sections_trait;
  ecs_int_t     nbr_zones_trait;

  ecs_int_t    *renum_sections;
  ecs_int_t    *renum_zones;
  char        **ptr_noms_sections;
  char        **ptr_noms_zones;

  bool         trait_section;
  bool         trait_zone;

  ecs_int_t    cpt_section;
  ecs_int_t    ind_section;
  ecs_int_t    ind_pass;
  ecs_int_t    ind_sub;
  ecs_int_t    ind_zone;
  ecs_int_t    ind_zs;

  char        *nom_grp = NULL;

  ecs_int_t   *tab_ind_zs_ent;
  ecs_int_t   *renum_zs;
  char       **ptr_nom_zs;

  /* Déclarations des variables de stockage        */
  /* avant transfert dans la structure du maillage */
  /*-----------------------------------------------*/

  bool         bool_aff_grp;
  bool         bool_termine;

  ecs_descr_t  * descr_grp;
  size_t         ent_cpt_elt[ECS_N_ENTMAIL];
  ecs_int_t    * ent_val_grp[ECS_N_ENTMAIL];

  ecs_table_t  * table_grp;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Initialisations */

  cel_dim = base_maillage->dim_entite;

  if (cel_dim == 1)
    return;

  else if (cel_dim == 2)
    ient_max = ECS_ENTMAIL_FAC;
  else
    ient_max = ECS_ENTMAIL_CEL;

  nbr_sections_trait = 0;
  nbr_zones_trait    = 0;

  renum_sections = NULL;
  renum_zones    = NULL;
  ptr_noms_sections = NULL;
  ptr_noms_zones    = NULL;

  trait_section = false;
  if (   tab_ind_section_cel != NULL
      || tab_ind_section_fac != NULL)
    trait_section = true;

  trait_zone = false;
  if (   tab_ind_zone_cel != NULL
      || tab_ind_zone_fac != NULL)
    trait_zone = true;


  if (trait_section == false && trait_zone == false)
    return;

  /* Renumérotation des noms et numéros de sections si nécessaire */

  if (trait_section == true) {

    cpt_section = 0;

    for (ind_zone = 0; ind_zone < nzones; ind_zone++) {
      ptr_zone = tab_zone + ind_zone;
      if (ptr_zone->type == CS_CG_ENUM(Unstructured))
        cpt_section += ptr_zone->nbr_sections;
    }

    ECS_MALLOC(renum_sections, cpt_section, ecs_int_t);
    ECS_MALLOC(ptr_noms_sections, cpt_section, char *);

    cpt_section = 0;

    for (ind_zone = 0; ind_zone < nzones; ind_zone++) {

      ptr_zone = tab_zone + ind_zone;

      if (ptr_zone->type == CS_CG_ENUM(Unstructured)) {

        for (ind_section = 0;
             ind_section < ptr_zone->nbr_sections;
             ind_section++) {

          ptr_section = ptr_zone->tab_sections + ind_section;

          for (ind_sub = 0; ind_sub < cpt_section; ind_sub++) {
            if (strcmp(ptr_noms_sections[ind_sub], ptr_section->nom) == 0)
              break;
          }

          nbr_sections_trait = ECS_MAX(ind_sub, nbr_sections_trait);

          renum_sections[cpt_section] = ind_sub;
          ptr_noms_sections[renum_sections[cpt_section]] = ptr_section->nom;

          cpt_section += 1;

        }

      }

    }

    nbr_sections_trait += 1;

  }

  /* Renumérotation des noms et numéros de zones si nécessaire */

  if (trait_zone == true) {

    ECS_MALLOC(renum_zones, nzones, ecs_int_t);
    ECS_MALLOC(ptr_noms_zones, nzones, char *);

    for (ind_zone = 0; ind_zone < nzones; ind_zone++) {

      ptr_zone = tab_zone + ind_zone;

      for (ind_sub = 0; ind_sub < ind_zone; ind_sub++) {
        if (strcmp(ptr_noms_zones[ind_sub], ptr_zone->nom) == 0)
          break;
      }

      nbr_zones_trait = ECS_MAX(ind_sub, nbr_zones_trait);

      renum_zones[ind_zone] = ind_sub;
      ptr_noms_zones[renum_zones[ind_zone]] = ptr_zone->nom;

    }

    nbr_zones_trait += 1;
  }

  /* Initialisations des tableaux locaux */

  for (ient = 0; ient < ECS_N_ENTMAIL; ient++) {
    ent_cpt_elt[ient] = 0;
    ent_val_grp[ient] = NULL;
  }

  /*-------------------------------------------------------*/
  /* Génération des groupes associés aux zones et sections */
  /*-------------------------------------------------------*/

  for (ind_pass = 0; ind_pass < 4; ind_pass++) {

    switch(ind_pass) {

    case 0:
      tab_ind_zs_ent = tab_ind_section_fac;
      ient = ient_max - 1;
      renum_zs = renum_sections;
      ptr_nom_zs = ptr_noms_sections;

      break;

    case 1:
      tab_ind_zs_ent = tab_ind_zone_fac;
      ient = ient_max - 1;
      renum_zs = renum_zones;
      ptr_nom_zs = ptr_noms_zones;

      break;

    case 2:
      tab_ind_zs_ent = tab_ind_section_cel;
      ient = ient_max;
      renum_zs = renum_sections;
      ptr_nom_zs = ptr_noms_sections;

      break;

    case 3:
      tab_ind_zs_ent = tab_ind_zone_cel;
      ient = ient_max;
      renum_zs = renum_zones;
      ptr_nom_zs = ptr_noms_zones;

      break;

    default:
      assert(0);

    }

    if (tab_ind_zs_ent != NULL && renum_zs != NULL && ient > ECS_ENTMAIL_NONE) {

      ind_zs = 0;

      do { /* Boucle tant que sur les groupes à créer */

        bool_termine = true;

        /* On alloue et initialise pour le groupe à lire */

        ent_cpt_elt[ient] = 0;

        if (ent_nbr_elt[ient] != 0)
          ECS_MALLOC(ent_val_grp[ient], ent_nbr_elt[ient], ecs_int_t);

        for (ind = 0; ind < (ecs_int_t)(ent_nbr_elt[ient]); ind++)
          ent_val_grp[ient][ind] = 0;

        for (ind = 0; ind < (ecs_int_t)(ent_nbr_elt[ient]); ind++) {

          if (renum_zs[tab_ind_zs_ent[ind]] == ind_zs) {

            ent_val_grp[ient][ind] = 1;
            ent_cpt_elt[ient]++;
            nom_grp = ptr_nom_zs[ind_zs];

          }
          else if (renum_zs[tab_ind_zs_ent[ind]] > ind_zs)

            bool_termine = false;

        }

        ind_zs += 1;

        /* Remplissage des entités du maillage */
        /*-------------------------------------*/

        bool_aff_grp = false;

        if (ent_cpt_elt[ient] != 0) {

          bool_aff_grp = true;

          assert(ent_cpt_elt[ient] <= (size_t)(ent_nbr_elt[ient]));

          /* Création du descripteur de table correspondant au groupe lu */

          descr_grp = ecs_descr__cree(ECS_DESCR_IDE_NUL,
                                      nom_grp);


          /* Transformation du tableau référencant le groupe en une table */

          table_grp = ecs_table__transforme_tableau(ent_nbr_elt[ient],
                                                    ent_val_grp[ient],
                                                    descr_grp);

          if (ent_table_grp[ient] != NULL) {
            ecs_table_att__assemble(ent_table_grp[ient],
                                    table_grp);
          }

          else
            ent_table_grp[ient] = table_grp;

        } /* Fin si le nombre d'éléments référencant le groupe n'est pas nul */


        /* Libération mémoire */

        if (ent_nbr_elt[ient] != 0)
          ECS_FREE(ent_val_grp[ient]);


        /* Affichage du bilan des données lues pour les groupes */
        /*------------------------------------------------------*/

        if (bool_aff_grp == true)
          printf("  %s \"%s\"\n",
                 _("Group"), nom_grp);

        ecs_maillage_pre__aff_nbr_par_ent(0,
                                          ent_cpt_elt,
                                          0);

        ent_cpt_elt[ient] = 0;


      } while (bool_termine == false);

    } /* Fin test si combinaison entité/section ou zones à traiter */

  } /* Fin boucle sur les combinaisons */

  /* Libération mémoire*/

  ECS_FREE(renum_sections);
  ECS_FREE(renum_zones);
  ECS_FREE(ptr_noms_sections);
  ECS_FREE(ptr_noms_zones);
}

/*----------------------------------------------------------------------------
 * Création des groupes basés sur les C.L.
 * et de manière optionnelle sur les zones et les sections
 *----------------------------------------------------------------------------*/

static void
ecs_loc_pre_cgns__cree_groupes(ecs_maillage_t              *maillage,
                               const ecs_loc_cgns_base_t   *base_maillage,
                               int                          nzones,
                               ecs_loc_cgns_zone_t         *tab_zone,
                               int                          nbr_boco_tot,
                               char                       **tab_nom_boco,
                               ecs_loc_cgns_boco_t         *tab_boco,
                               ecs_int_t                    nbr_boco_som,
                               ecs_int_t                  **ind_nom_boco_som,
                               ecs_int_t                  **nbr_sselt_boco_som,
                               ecs_int_t                 ***tab_sselt_boco_som,
                               ecs_int_t                   *tab_ind_section_cel,
                               ecs_int_t                   *tab_ind_section_fac,
                               ecs_int_t                   *tab_ind_zone_cel,
                               ecs_int_t                   *tab_ind_zone_fac)
{
  ecs_int_t     ient;

  size_t        ent_nbr_elt[ECS_N_ENTMAIL];

  ecs_table_t  *ent_table_grp[ECS_N_ENTMAIL];

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Initialisations */
  /*=================*/

  for (ient = 0; ient < ECS_N_ENTMAIL; ient++) {
    if (maillage->table_def[ient] != NULL)
      ent_nbr_elt[ient] = ecs_table__ret_elt_nbr(maillage->table_def[ient]);
    else
      ent_nbr_elt[ient] = 0;
    ent_table_grp[ient] = NULL;
  }

  /* Création des groupes */
  /*======================*/

  ecs_loc_pre_cgns__cree_grps_boco(base_maillage,
                                   tab_zone,
                                   nbr_boco_tot,
                                   tab_nom_boco,
                                   tab_boco,
                                   nbr_boco_som,
                                   ind_nom_boco_som,
                                   nbr_sselt_boco_som,
                                   tab_sselt_boco_som,
                                   ent_nbr_elt,
                                   ent_table_grp);

  ecs_loc_pre_cgns__cree_grps_zs(base_maillage,
                                 nzones,
                                 tab_zone,
                                 tab_ind_section_cel,
                                 tab_ind_section_fac,
                                 tab_ind_zone_cel,
                                 tab_ind_zone_fac,
                                 ent_nbr_elt,
                                 ent_table_grp);

  /* Transfert des valeurs lues dans les structures d'entité de maillage */
  /*=====================================================================*/

  for (ient = 0; ient < ECS_N_ENTMAIL; ient++) {
    assert(maillage->table_att[ient] == NULL);
    maillage->table_att[ient] = ent_table_grp[ient];
  }
}

/*============================================================================
 * Fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Lecture d'un fichier au format CGNS
 * et affectation des données dans la structure de maillage
 *----------------------------------------------------------------------------*/

ecs_maillage_t *
ecs_pre_cgns__lit_maillage(const char   *nom_fic_maillage,
                           int           num_maillage,
                           bool          cree_grp_cel_section,
                           bool          cree_grp_cel_zone,
                           bool          cree_grp_fac_section,
                           bool          cree_grp_fac_zone)
{
  int               ind;
  int               nbr_nom_boco;
  int               nbr_boco_tot;

  ecs_int_t         nbr_boco_som = 0;
  ecs_int_t        *ind_nom_boco_som;
  ecs_int_t        *nbr_sselt_boco_som;
  ecs_int_t       **tab_sselt_boco_som;

  char            **tab_nom_boco;

  ecs_int_t        *tab_ind_section_fac;
  ecs_int_t        *tab_ind_section_cel;
  ecs_int_t        *tab_ind_zone_cel;
  ecs_int_t        *tab_ind_zone_fac;

  ecs_int_t       **ptr_tab_ind_section_fac;
  ecs_int_t       **ptr_tab_ind_section_cel;
  ecs_int_t       **ptr_tab_ind_zone_cel;
  ecs_int_t       **ptr_tab_ind_zone_fac;

  ecs_loc_cgns_base_t  *base_maillage;
  ecs_loc_cgns_zone_t  *tab_zone;
  ecs_loc_cgns_boco_t  *tab_boco;

  /* Déclarations des variables pour CGNS */
  /*--------------------------------------*/

  char   nom_tmp[ECS_CGNS_TAILLE_NOM + 1];

  int    num_fic, num_base, cell_dim, phys_dim, nzones;
  int    ret_cgns = 0;
  float  version_cgns;

  bool ignore_vertex_bocos = false;
  if (cree_grp_fac_section || cree_grp_fac_zone)
    ignore_vertex_bocos = true;

  /* Création d'un maillage initialement vide (valeur de retour) */

  ecs_maillage_t  *maillage = ecs_maillage__cree_nodal();

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  printf(_("\n\n"
           "Reading mesh from file in CGNS format\n"
           "----------------------\n"));

  printf(_("  Mesh file: %s\n\n\n"),
         nom_fic_maillage);

  /* Ouverture du fichier CGNS en lecture */
  /*--------------------------------------*/

  base_maillage = ecs_loc_pre_cgns__cree(nom_fic_maillage,
                                         num_maillage);

  num_fic  = base_maillage->num_fic;
  num_base = base_maillage->num_base;

  ret_cgns = cg_version(num_fic, &version_cgns);

  if (ret_cgns == 0)
    printf(_("  CGNS version       : %3.1f\n"), version_cgns);

  /* Récupération titre si possible */

  ecs_loc_pre_cgns__aff_titre_cas(num_fic,
                                  num_base,
                                  _("  Title              : "));

  /* Autres informations sur la base */

  if (cg_base_read(num_fic, num_base, nom_tmp, &cell_dim, &phys_dim) != 0)
    ecs_error(__FILE__, __LINE__, 0,
              _("CGNS error:\n%s"), cg_get_error());

  printf(_("  Base name          : %s\n"), nom_tmp);
  printf(_("  Physical dimension : %d\n"), phys_dim);
  printf(_("  Cell dimension     : %d\n"), cell_dim);

  base_maillage->dim_entite = cell_dim;
  base_maillage->dim_espace = phys_dim;

  /* Nombre de zones */

  if (cg_nzones(num_fic, num_base, &nzones) != 0)
    ecs_error(__FILE__, __LINE__, 0,
              _("CGNS error:\n%s"), cg_get_error());
  printf(_("  Number of zones    : %d\n"), nzones);

  if (cell_dim < 2)
    ecs_error(__FILE__, __LINE__, 0,
              _("Mesh entities are of dimension: %d\n"),
              cell_dim);

  /* Lecture des informations principales sur les zones */
  /*----------------------------------------------------*/

  tab_zone = ecs_loc_pre_cgns__lit_zones(base_maillage,
                                         nzones);

  /* Lecture des conditions aux limites */
  /*------------------------------------*/

  ecs_loc_pre_cgns__lit_boco(base_maillage,
                             ignore_vertex_bocos,
                             nzones,
                             tab_zone,
                             &nbr_nom_boco,
                             &nbr_boco_tot,
                             &tab_nom_boco,
                             &tab_boco);

  /* Lecture des sommets */
  /*---------------------*/

  ecs_loc_pre_cgns__lit_som(maillage,
                            base_maillage,
                            nzones,
                            tab_zone);

  /* Lecture des éléments */
  /*----------------------*/

  tab_ind_section_fac = NULL;
  tab_ind_section_cel = NULL;
  tab_ind_zone_cel    = NULL;
  tab_ind_zone_fac    = NULL;

  ptr_tab_ind_section_cel
    = (cree_grp_cel_section == true) ? &tab_ind_section_cel : NULL;
  ptr_tab_ind_section_fac
    = (cree_grp_fac_section == true) ? &tab_ind_section_fac : NULL;
  ptr_tab_ind_zone_cel
    = (cree_grp_cel_zone == true) ? &tab_ind_zone_cel : NULL;
  ptr_tab_ind_zone_fac
    = (cree_grp_fac_zone == true) ? &tab_ind_zone_fac : NULL;

  ecs_loc_pre_cgns__lit_ele(maillage,
                            base_maillage,
                            nzones,
                            tab_zone,
                            nbr_boco_tot,
                            tab_boco,
                            &nbr_boco_som,
                            &ind_nom_boco_som,
                            &nbr_sselt_boco_som,
                            &tab_sselt_boco_som,
                            ptr_tab_ind_section_cel,
                            ptr_tab_ind_section_fac,
                            ptr_tab_ind_zone_cel,
                            ptr_tab_ind_zone_fac);

  /* Création des groupes */
  /*----------------------*/

  ecs_loc_pre_cgns__cree_groupes(maillage,
                                 base_maillage,
                                 nzones,
                                 tab_zone,
                                 nbr_boco_tot,
                                 tab_nom_boco,
                                 tab_boco,
                                 nbr_boco_som,
                                 &ind_nom_boco_som,
                                 &nbr_sselt_boco_som,
                                 &tab_sselt_boco_som,
                                 tab_ind_section_cel,
                                 tab_ind_section_fac,
                                 tab_ind_zone_cel,
                                 tab_ind_zone_fac);


  /* Libération mémoire */

  ECS_FREE(tab_ind_section_cel);
  ECS_FREE(tab_ind_section_fac);
  ECS_FREE(tab_ind_zone_cel);
  ECS_FREE(tab_ind_zone_fac);

  for (ind = 0; ind < nbr_nom_boco; ind++)
    ECS_FREE(tab_nom_boco[ind]);
  ECS_FREE(tab_nom_boco);

  ECS_FREE(tab_boco);

  for (ind = 0; ind < nzones; ind++) {
    if ((tab_zone + ind)->tab_sections != NULL)
      ECS_FREE((tab_zone + ind)->tab_sections);
    if ((tab_zone + ind)->renum != NULL)
      ECS_FREE((tab_zone + ind)->renum);
  }
  ECS_FREE(tab_zone);

  /* Fermeture du fichier de lecture du maillage */

  ecs_loc_pre_cgns__detruit(base_maillage);

  /* Renvoi de la structure de maillage */

  return maillage;
}

/*----------------------------------------------------------------------------*/

#endif /* HAVE_CGNS */
