/*============================================================================
 *  Définitions des fonctions
 *   associées à la structure `ecs_champ_t' décrivant un champ
 *   et réalisant les sorties pour post-traitement CGNS
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

#if defined(HAVE_CGNS)

/*============================================================================*
 *                                 Visibilité
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie standard C
 *---------------------------------------------------------------------------*/

#include <assert.h>
#include <string.h>

#include <cgnslib.h>


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage global "Utilitaire"
 *---------------------------------------------------------------------------*/

#include "ecs_def.h"
#include "ecs_elt_typ_liste.h"
#include "ecs_fic.h"
#include "ecs_mem.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage global "Post-Traitement"
 *---------------------------------------------------------------------------*/

#include "ecs_post.h"
#include "ecs_post_cgns.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles des paquetages visibles
 *---------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage courant
 *---------------------------------------------------------------------------*/

#include "ecs_champ.h"


/*----------------------------------------------------------------------------
 *  Fichier  `include' du  paquetage courant associé au fichier courant
 *---------------------------------------------------------------------------*/

#include "ecs_post_cgns.h"
#include "ecs_champ_post_cgns.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' privés   du  paquetage courant
 *---------------------------------------------------------------------------*/

#include "ecs_post_cgns_priv.h"
#include "ecs_champ_priv.h"


/*============================================================================
 * Définitions de paramètres et macros
 *============================================================================*/

/* Compatibilité avec diverses versions de CGNS */

#if defined(CGNS_SCOPE_ENUMS)
#define CS_CG_ENUM(e) CG_ ## e
#else
#define CS_CG_ENUM(e) e
#endif


/*============================================================================
 *                              Fonctions privées
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction qui renvoie un pointeur sur une sous-structure associée
 *  à un maillage pour un cas de sortie CGNS.
 *---------------------------------------------------------------------------*/

static ecs_post_cgns_base_t *
ecs_loc_champ_post_cgns__base(const ecs_post_cgns_t  *cas_cgns,
                              const char             *nom_maillage)
{
  ecs_int_t  ind;
  ecs_post_cgns_base_t  *base_cgns = NULL;

  /* Recherche du maillage */
  /*-----------------------*/

  for (ind = 0; ind < cas_cgns->nbr_bases; ind++) {
    base_cgns = cas_cgns->tab_bases[ind];
    if (strcmp(nom_maillage, base_cgns->nom_maillage) == 0)
      break;
  }

  if (ind >= cas_cgns->nbr_bases)
    ecs_error(__FILE__, __LINE__, 0,
              _("CGNS: no mesh named \"%s\".\n"
                "is associated to file: \"%s\"\n"),
              nom_maillage, base_cgns->nom_fic);

  /* Réouverture du fichier associé si nécessaire */

  if (base_cgns->fic_ouvert == false) {

    if (   cg_open(base_cgns->nom_fic, CG_MODE_MODIFY, &(base_cgns->num_fic))
        != CG_OK)
      ecs_error(__FILE__, __LINE__, 0,
                _("CGNS: error re-opening file \"%s\":\n%s"),
                base_cgns->nom_fic, cg_get_error());

    base_cgns->fic_ouvert = true;

  }
  return  base_cgns;
}

/*----------------------------------------------------------------------------
 * Write per-element values for a given array.
 *
 * Values corresponding to polyhedra are ignored.
 *---------------------------------------------------------------------------*/

static void
_champ_post_cgns__ecr_fam(const int        *elt_fam,
                          const char       *nom_maillage,
                          ecs_post_cgns_t  *cas_cgns)
{
  ecs_post_cgns_base_t  *base_cgns;

  /* Declarations des variables pour CGNS */
  /*-------------------------------------*/

  char  nom_champ_cgns[ECS_CGNS_TAILLE_NOM + 1];

  CS_CG_ENUM(DataType_t)  type_val_cgns;

  int  num_sol;
  int  num_champ;
  int  ret_cgns;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(elt_fam != NULL);

  /* Recherche de la base CGNS */
  /*---------------------------*/

  base_cgns = ecs_loc_champ_post_cgns__base(cas_cgns,
                                            nom_maillage);

  if (base_cgns == NULL)
    return;

  /* Ajout d'une solution au cas CGNS si nécessaire */
  /*------------------------------------------------*/

  if (base_cgns->num_sol_stat < 0) {

    if (cg_sol_write(base_cgns->num_fic,
                     1,
                     1,
                     _("Mesh_values"),
                     CS_CG_ENUM(CellCenter),
                     &num_sol) != CG_OK)

      ecs_error(__FILE__, __LINE__, 0,
                _("CGNS: error writing a solution node\n"
                  "File name: \"%s\"; base : \"1\"\n%s"),
                base_cgns->nom_fic, cg_get_error());

    base_cgns->num_sol_stat = num_sol;
  }

  /* Nom du champ */

  strncpy(nom_champ_cgns, _("Family"), ECS_CGNS_TAILLE_NOM);
  nom_champ_cgns[ECS_CGNS_TAILLE_NOM] = '\0';

  /* Écriture du champ */
  /*-------------------*/

  type_val_cgns = CS_CG_ENUM(Integer);

  ret_cgns = cg_field_write(base_cgns->num_fic,
                            1,
                            1,
                            base_cgns->num_sol_stat,
                            type_val_cgns,
                            nom_champ_cgns,
                            elt_fam,
                            &num_champ);

  if (ret_cgns != CG_OK)
    ecs_error(__FILE__, __LINE__, 0,
              _("CGNS: error writing field \"%s\"\n%s"),
              nom_champ_cgns, cg_get_error());
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Write element connectivity based on geometric type.
 *
 * Elements must have been previously sorted by type, and polyhedra
 * are ignored.
 *---------------------------------------------------------------------------*/

void
ecs_champ_post_cgns__ecr_connect(const char            *nom_maillage,
                                 size_t                 n_vertices,
                                 const ecs_coord_t      vertex_coords[],
                                 ecs_champ_t           *champ_def,
                                 const int              elt_fam[],
                                 const ecs_tab_int_t   *tab_elt_typ_geo,
                                 ecs_post_cgns_t       *cas_cgns)
{
  size_t     ielt;
  size_t     ival_deb;
  size_t     isom;
  int        icoo;
  ecs_int_t  ind;
  ecs_int_t  ind_typ;

  size_t     cpt_elt;
  size_t     cpt_elt_fin;
  size_t     nbr_elt;
  size_t     nbr_som;
  size_t     nbr_som_elt;
  size_t     nbr_val;
  int        elt_typ_ref;

  ecs_size_t * def_pos_tab;
  ecs_int_t  * def_val_tab;

  int         cpt_section;
  size_t      nbr_elt_typ;
  int         isize[3];
  int         num_coord;
  int         num_section;
  int         num_zone;
  int         ret_cgns;

  int         type_cgns[ECS_ELT_TYP_FIN];
  int         type_cgns_loc;

  int        *def_elt;

  double     *coo_temp;

  char        nom_section[32 + 1];

  char const  *nom_coord[3] = {"CoordinateX",
                               "CoordinateY",
                               "CoordinateZ"};

  ecs_post_cgns_base_t  *base_cgns;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(vertex_coords != NULL);
  assert(champ_def != NULL);

  assert(cas_cgns != NULL);

  /* Recherche de la base CGNS */
  /*---------------------------*/

  base_cgns = ecs_loc_champ_post_cgns__base(cas_cgns,
                                            nom_maillage);

  /* Dimensions */
  /*------------*/

  nbr_som = n_vertices;
  nbr_elt = champ_def->nbr;

  isize[0] = nbr_som;
  isize[1] = nbr_elt - tab_elt_typ_geo->val[ECS_ELT_TYP_CEL_POLY];

  isize[2] = 0;       /* éléments de bord non triés */

  ret_cgns = cg_zone_write(base_cgns->num_fic,
                           1,
                           "Zone 1",
                           isize,
                           CS_CG_ENUM(Unstructured),
                           &num_zone);

  if (ret_cgns != CG_OK)
    ecs_error(__FILE__, __LINE__, 0,
              _("CGNS: error writing file \"%s\":\n"
                "Name of mesh to write: \"%s\"\n%s"),
              base_cgns->nom_fic, base_cgns->nom_maillage, cg_get_error());

  assert (num_zone == 1);

  /* Écriture des coordonnées */
  /*--------------------------*/

  /* Remplissage du tableau temporaire et écriture */

  ECS_MALLOC(coo_temp, nbr_som, double);

  for (icoo = 0; icoo < base_cgns->dim_espace; icoo++) {

    for (isom = 0; isom < nbr_som; isom++)
      coo_temp[isom] = vertex_coords[isom * 3 + icoo];

    ret_cgns = cg_coord_write(base_cgns->num_fic,
                              1,
                              1,
                              CS_CG_ENUM(RealDouble),
                              nom_coord[icoo],
                              coo_temp,
                              &num_coord);

    if (ret_cgns != CG_OK)
      ecs_error(__FILE__, __LINE__, 0,
                _("CGNS: error writing coordinates\n"
                  "Name of mesh to write: \"%s\"\n%s"),
                base_cgns->nom_maillage, cg_get_error());

  }

  ECS_FREE(coo_temp);

  /* Écriture des éléments */
  /*------------------------*/

  ecs_champ__regle_en_pos(champ_def);

  def_pos_tab = champ_def->pos;
  def_val_tab = champ_def->val;

  cpt_section = 0;

  cpt_elt = 0;

  elt_typ_ref = -1;

  for (ind_typ = ECS_ELT_TYP_NUL; ind_typ < ECS_ELT_TYP_FIN; ind_typ++)
    type_cgns[ind_typ] = -2;

  type_cgns[ECS_ELT_TYP_NUL] = -1;
  type_cgns[ECS_ELT_TYP_FAC_TRIA] = CS_CG_ENUM(TRI_3);
  type_cgns[ECS_ELT_TYP_FAC_QUAD] = CS_CG_ENUM(QUAD_4);
  type_cgns[ECS_ELT_TYP_CEL_TETRA] = CS_CG_ENUM(TETRA_4);
  type_cgns[ECS_ELT_TYP_CEL_PYRAM] = CS_CG_ENUM(PYRA_5);
  type_cgns[ECS_ELT_TYP_CEL_PRISM] = CS_CG_ENUM(PENTA_6);
  type_cgns[ECS_ELT_TYP_CEL_HEXA] = CS_CG_ENUM(HEXA_8);
  type_cgns[ECS_ELT_TYP_FAC_POLY] = CS_CG_ENUM(NGON_n);
  type_cgns[ECS_ELT_TYP_CEL_POLY] = - 1;

#if defined(DEBUG) && !defined(NDEBUG)
  for (ind_typ = ECS_ELT_TYP_NUL; ind_typ < ECS_ELT_TYP_FIN; ind_typ++)
    assert(type_cgns[ind_typ] != -2);
#endif

  while (cpt_elt < nbr_elt) {

    /* Recherche du prochain type d'élément utilisé */

    elt_typ_ref += 1;

    while (tab_elt_typ_geo->val[elt_typ_ref] == 0)
      elt_typ_ref++;

    nbr_elt_typ = tab_elt_typ_geo->val[elt_typ_ref];
    cpt_elt_fin = cpt_elt + nbr_elt_typ;

    /* Création des connectivités */
    /*----------------------------*/

    ind = 0;
    nbr_val = def_pos_tab[cpt_elt_fin] - def_pos_tab[cpt_elt];

    if (   elt_typ_ref != ECS_ELT_TYP_FAC_POLY
        && elt_typ_ref != ECS_ELT_TYP_CEL_POLY) {

      type_cgns_loc = type_cgns[elt_typ_ref];

      ECS_MALLOC(def_elt, nbr_val, int);

      for (ielt = cpt_elt; ielt < cpt_elt_fin; ielt++) {

        ival_deb = def_pos_tab[ielt] - 1;

        nbr_som_elt = def_pos_tab[ielt + 1] - def_pos_tab[ielt];

        /* Numérotation locale CGNS identique à la numéoration
           locale interne pour les éléments linéaires */

        for (isom = 0; isom < nbr_som_elt; isom++)
          def_elt[ind++] = def_val_tab[ival_deb + isom];

      }

    }
    else if (elt_typ_ref == ECS_ELT_TYP_FAC_POLY) { /* Cas des polygones */

#if 0

      type_cgns_loc = CS_CG_ENUM(MIXED);

      ECS_MALLOC(def_elt, nbr_val + cpt_elt_fin - cpt_elt, int);

      for (ielt = cpt_elt; ielt < cpt_elt_fin; ielt++) {

        def_elt[ind++]
          = def_pos_tab[ielt + 1] - def_pos_tab[ielt] + CS_CG_ENUM(NGON_n);

        for (ival = def_pos_tab[ielt    ] - 1;
             ival < def_pos_tab[ielt + 1] - 1;
             ival++)
          def_elt[ind++] = def_val_tab[ival];

      }

#else

      cpt_elt += nbr_elt_typ;

      ecs_warn();
      printf(_("CGNS: in mesh: \"%s\",\n"
               "%d polygonal faces are ignored.\n"),
             base_cgns->nom_maillage, (int)nbr_elt_typ);

      break;

#endif

    }

    else { /* Cas des polyèdres (ignorés) */

      cpt_elt += nbr_elt_typ;

      ecs_warn();
      printf(_("CGNS: in mesh: \"%s\",\n"
               "%d polyhedral cells are ignored.\n"),
             base_cgns->nom_maillage, (int)nbr_elt_typ);

      break;

    }

    /* Écriture de la section (sections séparées) */

    cpt_section += 1;
    sprintf(nom_section, "Section %2d", cpt_section);

    ret_cgns = cg_section_write(base_cgns->num_fic,
                                1,
                                1,
                                nom_section,
                                type_cgns_loc,
                                (int)(cpt_elt + 1),
                                (int)(cpt_elt + nbr_elt_typ),
                                0,
                                def_elt,
                                &num_section);

    ECS_FREE(def_elt);

    /* On s'apprête à passer au type d'élément suivant */

    cpt_elt += nbr_elt_typ;
  }

  /* Nettoyage avant la sortie */
  /*---------------------------*/

  ecs_champ__libere_pos_tab(champ_def, def_pos_tab);

  /* Familles des éléments */
  /*-----------------------*/

  if (elt_fam != NULL)
    _champ_post_cgns__ecr_fam(elt_fam,
                              nom_maillage,
                              cas_cgns);
}

/*----------------------------------------------------------------------------
 * Write per-element values for a given array.
 *
 * Values corresponding to polyhedra are ignored.
 *---------------------------------------------------------------------------*/

void
ecs_champ_post_cgns__ecr_val(const ecs_tab_int_t  *tab_val,
                             const char           *nom_maillage,
                             const char           *nom_champ,
                             ecs_post_cgns_t      *cas_cgns)
{
  ecs_post_cgns_base_t  *base_cgns;

  /* Declarations des variables pour CGNS */
  /*-------------------------------------*/

  char               nom_champ_cgns[ECS_CGNS_TAILLE_NOM + 1];

  CS_CG_ENUM(DataType_t)  type_val_cgns;

  int                num_sol;
  int                num_champ;
  int                ret_cgns;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(tab_val != NULL);

  /* Recherche de la base CGNS */
  /*---------------------------*/

  base_cgns = ecs_loc_champ_post_cgns__base(cas_cgns,
                                            nom_maillage);

  if (base_cgns == NULL)
    return;

  /* Ajout d'une solution au cas CGNS si nécessaire */
  /*------------------------------------------------*/

  if (base_cgns->num_sol_stat < 0) {

    if (cg_sol_write(base_cgns->num_fic,
                     1,
                     1,
                     _("Mesh_values"),
                     CS_CG_ENUM(CellCenter),
                     &num_sol) != CG_OK)

      ecs_error(__FILE__, __LINE__, 0,
                _("CGNS: error writing a solution node\n"
                  "File name: \"%s\"; base : \"1\"\n%s"),
                base_cgns->nom_fic, cg_get_error());

    base_cgns->num_sol_stat = num_sol;

  }

  /* Nom du champ */

  strncpy(nom_champ_cgns, nom_champ, ECS_CGNS_TAILLE_NOM);
  nom_champ_cgns[ECS_CGNS_TAILLE_NOM] = '\0';

  /* Écriture du champ */
  /*-------------------*/

  assert(sizeof(ecs_int_t) == sizeof(int));

  type_val_cgns = CS_CG_ENUM(Integer);

  if (tab_val->nbr > 0) {

    ret_cgns = cg_field_write(base_cgns->num_fic,
                              1,
                              1,
                              base_cgns->num_sol_stat,
                              type_val_cgns,
                              nom_champ_cgns,
                              (void *)tab_val->val,
                              &num_champ);

    if (ret_cgns != CG_OK)
      ecs_error(__FILE__, __LINE__, 0,
                _("CGNS: error writing field \"%s\"\n%s"),
                nom_champ_cgns, cg_get_error());

  }
}

#endif /* HAVE_CGNS */

/*----------------------------------------------------------------------------*/

