/*============================================================================
 *  Définition de la fonction
 *   de lecture d'un fichier de maillage au format Gmsh
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie standard C
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>


/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage global "Utilitaire"
 *----------------------------------------------------------------------------*/

#include "ecs_def.h"
#include "ecs_elt_typ_liste.h"
#include "ecs_file.h"
#include "ecs_mem.h"
#include "ecs_tab.h"


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
 *  Fichier  `include' du  paquetage courant associé au fichier courant
 *----------------------------------------------------------------------------*/

#include "ecs_pre_gmsh.h"


/*----------------------------------------------------------------------------
 *  Fichiers `include' privés   du  paquetage courant
 *----------------------------------------------------------------------------*/


/*============================================================================
 *  Définitions de paramètres-macros
 *============================================================================*/

#define ECS_GMSH_NBR_TYP_ELT        31 /* Number of Gmsh element types */
#define ECS_GMSH_NBR_MAX_SOM        56 /* Max number of nodes per element */

/* Valeurs associées aux dimensionnements pour la lecture des lignes */

#define ECS_LOC_LNG_MAX_CHAINE_GMSH  2048 /* Max file line length */

typedef int gmsh_int_t;

/*============================================================================
 *                          Définitions de types
 *============================================================================*/

typedef enum {
  GMSH_SEG2 = 1, /*  1 */
  GMSH_TRIA3,    /*  2 */
  GMSH_QUAD4,    /*  3 */
  GMSH_TETRA4,   /*  4 */
  GMSH_HEXA8,    /*  5 */
  GMSH_PENTA6,   /*  6 */
  GMSH_PYRA5,    /*  7 */
  GMSH_SEG3,     /*  8 */
  GMSH_TRIA6,    /*  9 */
  GMSH_QUAD9,    /* 10 */
  GMSH_TETRA10,  /* 11 */
  GMSH_HEXA27,   /* 12 */
  GMSH_PENTA18,  /* 13 */
  GMSH_PYRA14,   /* 14 */
  GMSH_POINT1,   /* 15 */
  GMSH_QUAD8,    /* 16 */
  GMSH_HEXA20,   /* 17 */
  GMSH_PENTA15,  /* 18 */
  GMSH_PYRA13,   /* 19 */
  GMSH_TRIA9,    /* 20 */
  GMSH_TRIA10,   /* 21 */
  GMSH_TRIA12,   /* 22 */
  GMSH_TRIA15,   /* 23 */
  GMSH_TRIA15I,  /* 24 */
  GMSH_TRIA21,   /* 25 */
  GMSH_SEG4,     /* 26 */
  GMSH_SEG5,     /* 27 */
  GMSH_SEG6,     /* 28 */
  GMSH_TETRA20,  /* 29 */
  GMSH_TETRA35,  /* 30 */
  GMSH_TETRA56   /* 31 */
} ecs_gmsh_elt_typ_t;

typedef struct {

  ecs_gmsh_elt_typ_t     gmsh_typ  ; /* Type Gmsh de l'élément */
  ecs_elt_typ_t          ecs_typ   ; /* Type ECS de l'élément */
  ecs_int_t              nbr_som   ; /* Nombre de sommets */
                                      /* Liste des numéros de sommet ECS */
  ecs_int_t              num_som[ECS_GMSH_NBR_MAX_SOM];

} ecs_gmsh_elt_t;

/*============================================================================
 *  Définitions de variables globales statiques
 *============================================================================*/

const ecs_gmsh_elt_t ecs_gmsh_elt_liste_c[ECS_GMSH_NBR_TYP_ELT] = {
  {                                /* 1 */
    GMSH_SEG2,
    ECS_ELT_TYP_NUL,
    2,
    { 0 }
  },
  {                                /* 2 */
    GMSH_TRIA3,
    ECS_ELT_TYP_FAC_TRIA,
    3,
    { 1, 2, 3 }
  },
  {                                /* 3 */
    GMSH_QUAD4,
    ECS_ELT_TYP_FAC_QUAD,
    4,
    { 1, 2, 3, 4 }
  },
  {                                /* 4 */
    GMSH_TETRA4,
    ECS_ELT_TYP_CEL_TETRA,
    4,
    { 1, 2, 3, 4 }
  },
  {                                /* 5 */
    GMSH_HEXA8,
    ECS_ELT_TYP_CEL_HEXA,
    8,
    { 1, 2, 3, 4, 5, 6, 7, 8 },
  },
  {                                /* 6 */
    GMSH_PENTA6,
    ECS_ELT_TYP_CEL_PRISM,
    6,
    { 1, 2, 3, 4, 5, 6 }
  },
  {                                /* 7 */
    GMSH_PYRA5,
    ECS_ELT_TYP_CEL_PYRAM,
    5,
    { 1, 2, 3, 4, 5 }
  },
  {                                /* 8 */
    GMSH_SEG3,
    ECS_ELT_TYP_NUL,
    3,
    { 0 }
  },
  {                                /* 9 */
    GMSH_TRIA6,
    ECS_ELT_TYP_FAC_TRIA,
    6,
    { 1, 2, 3 }
  },
  {                               /* 10 */
    GMSH_QUAD9,
    ECS_ELT_TYP_FAC_QUAD,
    9,
    { 1, 2, 3, 4 }
  },
  {                               /* 11 */
    GMSH_TETRA10,
    ECS_ELT_TYP_CEL_TETRA,
    10,
    { 1, 2, 3, 4 }
  },
  {                               /* 12 */
    GMSH_HEXA27,
    ECS_ELT_TYP_CEL_HEXA,
    27,
    { 1, 2, 3, 4, 5, 6, 7, 8 },
  },
  {                               /* 13 */
    GMSH_PENTA18,
    ECS_ELT_TYP_CEL_PRISM,
    18,
    { 1, 2, 3, 4, 5, 6 }
  },
  {                               /* 14 */
    GMSH_PYRA14,
    ECS_ELT_TYP_CEL_PYRAM,
    14,
    { 1, 2, 3, 4, 5 }
  },
  {                               /* 15 */
    GMSH_POINT1,
    ECS_ELT_TYP_NUL,
    1,
    { 0 }
  },
  {                               /* 16 */
    GMSH_QUAD8,
    ECS_ELT_TYP_FAC_QUAD,
    8,
    { 1, 2, 3, 4 }
  },
  {                               /* 17 */
    GMSH_HEXA20,
    ECS_ELT_TYP_CEL_HEXA,
    20,
    { 1, 2, 3, 4, 5, 6, 7, 8 },
  },
  {                               /* 18 */
    GMSH_PENTA15,
    ECS_ELT_TYP_CEL_PRISM,
    15,
    { 1, 2, 3, 4, 5, 6 }
  },
  {                               /* 19 */
    GMSH_PYRA13,
    ECS_ELT_TYP_CEL_PYRAM,
    13,
    { 1, 2, 3, 4, 5 }
  },
  {                               /* 20 */
    GMSH_TRIA9,
    ECS_ELT_TYP_FAC_TRIA,
    9,
    { 1, 2, 3 }
  },
  {                               /* 21 */
    GMSH_TRIA10,
    ECS_ELT_TYP_FAC_TRIA,
    10,
    { 1, 2, 3 }
  },
  {                               /* 22 */
    GMSH_TRIA12,
    ECS_ELT_TYP_FAC_TRIA,
    12,
    { 1, 2, 3 }
  },
  {                               /* 23 */
    GMSH_TRIA15,
    ECS_ELT_TYP_FAC_TRIA,
    15,
    { 1, 2, 3 }
  },
  {                               /* 24 */
    GMSH_TRIA15I,
    ECS_ELT_TYP_FAC_TRIA,
    15,
    { 1, 2, 3 }
  },
  {                               /* 26 */
    GMSH_SEG4,
    ECS_ELT_TYP_NUL,
    4,
    { 0 }
  },
  {                               /* 27 */
    GMSH_SEG5,
    ECS_ELT_TYP_NUL,
    5,
    { 0 }
  },
  {                               /* 28 */
    GMSH_SEG6,
    ECS_ELT_TYP_NUL,
    6,
    { 0 }
  },
  {                               /* 29 */
    GMSH_TETRA20,
    ECS_ELT_TYP_CEL_TETRA,
    20,
    { 1, 2, 3, 4 }
  },
  {                               /* 30 */
    GMSH_TETRA35,
    ECS_ELT_TYP_CEL_TETRA,
    35,
    { 1, 2, 3, 4 }
  },
  {                               /* 31 */
    GMSH_TETRA56,
    ECS_ELT_TYP_CEL_TETRA,
    56,
    { 1, 2, 3, 4 }
  }
};

/*============================================================================
 *  Fonctions privées
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Lecture et vérification de la version du format (pour version 2.0)
 *----------------------------------------------------------------------------*/

static void
ecs_loc_pre_gmsh__lit_vers_format(ecs_file_t     *fic_maillage,
                                  int            *version,
                                  int            *type,
                                  int            *num_ligne)
{
  char  chaine[ECS_LOC_LNG_MAX_CHAINE_GMSH];
  int   retour;
  int   taille_coo;
  float fversion;

  /* Décodage de la chaine de version */

  ecs_file_gets(chaine, ECS_LOC_LNG_MAX_CHAINE_GMSH, fic_maillage, num_ligne);

  retour = sscanf(chaine,"%f %d %d", &fversion, type, &taille_coo);

  if (retour != 3)
    ecs_error(__FILE__, __LINE__, 0,
              _("The Gmsh version specification of file\n"
                "\"%s\" is invalid or unreadable.)"));

  if (*type != 0 && *type != 1)
    ecs_error(__FILE__, __LINE__, 0,
              _("The Gmsh type specification of file\n"
                "\"%s\" is neither 0 (text file) nor 1 (binary file).\n"
                "This case is currently not handled."));

  printf(_("  Gmsh format version:     %2.1f\n"
           "  Size given for real numbers: %d\n\n"),
         fversion, taille_coo);

  *version = (int) (fversion*10 + 1e-3);

  if (*type == 1) {

    int   un;

    ecs_file_set_type(fic_maillage, ECS_FILE_TYPE_BINARY);
    ecs_file_read(&un, sizeof(gmsh_int_t), 1, fic_maillage);
    if (un != 1)
      ecs_file_set_swap_endian(fic_maillage, 1);

    ecs_file_set_type(fic_maillage, ECS_FILE_TYPE_TEXT);
    ecs_file_gets(chaine, 2, fic_maillage, num_ligne);

  }

  /* Ligne de fin de rubrique */

  ecs_file_gets(chaine, ECS_LOC_LNG_MAX_CHAINE_GMSH, fic_maillage, num_ligne);

  if (strncmp(chaine, "$EndMeshFormat", strlen("$EndMeshFormat")) != 0)
    ecs_error(__FILE__, __LINE__, 0,
              _("The Gmsh format version specification end for file\n"
                "\"%s\" is not present at the expected place.)"));
}

/*----------------------------------------------------------------------------
 *  Lecture des coordonnées des noeuds
 *----------------------------------------------------------------------------*/

static void
ecs_loc_pre_gmsh__lit_nodes(ecs_maillage_t   *maillage,
                            ecs_file_t       *fic_maillage,
                            int               type_fmt_gmsh,
                            int              *num_ligne,
                            ecs_int_t       **som_val_label)
{
  char   chaine[ECS_LOC_LNG_MAX_CHAINE_GMSH];
  int    retour;

  /* Variables Gmsh lues */

  long         label;
  double       coord[3];
  ecs_int_t    nbr_nod;

  ecs_int_t    icoo;

  /* Stockage avant transfert */

  ecs_int_t    ind_nod = 0;

  ecs_coord_t * som_val_coord;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Initialisations */
  /*=================*/

  /* Décodage du nombre de noeuds */

  ecs_file_gets(chaine, ECS_LOC_LNG_MAX_CHAINE_GMSH, fic_maillage, num_ligne);

  nbr_nod = (ecs_int_t)(atol(chaine));

  /* Allocation mémoire associée */

  ECS_MALLOC(som_val_coord, nbr_nod * 3, ecs_coord_t);
  ECS_MALLOC(*som_val_label, nbr_nod, ecs_int_t);

  /* Boucle de lecture des noeuds */
  /*==============================*/

  if (type_fmt_gmsh == 0) {

    for (ind_nod = 0; ind_nod < nbr_nod; ind_nod++) {

      ecs_file_gets(chaine,
                    ECS_LOC_LNG_MAX_CHAINE_GMSH,
                    fic_maillage,
                    num_ligne);

      retour = sscanf(chaine,"%ld %lg %lg %lg",
                      &label, coord, coord+1, coord+2);

      if (retour != 4)
        ecs_error(__FILE__, __LINE__, 0,
                  _("Error decoding line %ld of file\n\"%s\" :\n"
                    "The description of point <%ld> was expected in the form "
                    "\"label x y z\"."),
                  (long)(*num_ligne), ecs_file_get_name(fic_maillage),
                  (long)(ind_nod+1));

      /* Étiquette du noeud lu */

      (*som_val_label)[ind_nod] = (ecs_int_t)label;

      /* Coordonnées du noeud lu */

      for (icoo = 0; icoo < 3; icoo++)
        som_val_coord[ind_nod * 3 + icoo] = (ecs_coord_t)(coord[icoo]);

    }

  }
  else {

    ecs_file_set_type(fic_maillage, ECS_FILE_TYPE_BINARY);

    for (ind_nod = 0; ind_nod < nbr_nod; ind_nod++) {

      gmsh_int_t _label;
      double xyz[3];

      ecs_file_read(&_label, sizeof(gmsh_int_t), 1, fic_maillage);
      ecs_file_read(&xyz, sizeof(double), 3, fic_maillage);

      (*som_val_label)[ind_nod] = (ecs_int_t)_label;
      for (icoo = 0; icoo < 3; icoo++)
        som_val_coord[ind_nod * 3 + icoo] = xyz[icoo];

    }

    ecs_file_set_type(fic_maillage, ECS_FILE_TYPE_TEXT);

    ecs_file_gets(chaine, 2, fic_maillage, num_ligne);

  }

  /* Ligne de fin de rubrique */

  ecs_file_gets(chaine, ECS_LOC_LNG_MAX_CHAINE_GMSH, fic_maillage, num_ligne);

  if (   (strncmp(chaine, "$EndNodes", strlen("$EndNodes")) != 0)
      && (strncmp(chaine, "$ENDNOD", strlen("$ENDNOD")) != 0))
    ecs_error(__FILE__, __LINE__, 0,
              _("The end of the points list specification of file\n"
                "\"%s\" is not present at the expected place.)"),
              ecs_file_get_name(fic_maillage));

  /* Transfert des valeurs lues dans la structure d'entité de maillage */
  /*===================================================================*/

  ecs_maillage_pre__cree_som(maillage,
                             nbr_nod,
                             som_val_coord);
}

/*----------------------------------------------------------------------------
 * Read mesh physical names
 *----------------------------------------------------------------------------*/

static void
_read_physical_names(ecs_file_t       *mesh_file,
                     int              *line_num,
                     int              *n_phys_names,
                     int              *n_phys_names_max,
                     int             **phys_name_id,
                     char           ***phys_name)
{
  char   line[ECS_LOC_LNG_MAX_CHAINE_GMSH];
  char   name[ECS_LOC_LNG_MAX_CHAINE_GMSH];
  int    retval;

  int    cpt_phys_names = 0;
  int    _n_phys_names_max = -1;

  int    _n_phys_names = 0;
  int   *_phys_name_id = NULL;
  char **_phys_name = NULL;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  ecs_file_gets(line, ECS_LOC_LNG_MAX_CHAINE_GMSH, mesh_file, line_num);

  retval = sscanf(line,"%d", &_n_phys_names);

  if (retval != 1)
    ecs_error(__FILE__, __LINE__, 0,
              _("Error decoding line %ld of file\n\"%s\" :\n"
                "The information on physical groups form "
                "\"numPhysicalNames\"."),
                (long)(*line_num), ecs_file_get_name(mesh_file));

  for (int i = 0; i < _n_phys_names; i++) {

    gmsh_int_t  dim_phys_name, tag_phys_name;

    ecs_file_gets(line,
                  ECS_LOC_LNG_MAX_CHAINE_GMSH,
                  mesh_file,
                  line_num);

    retval = sscanf(line, "%d %d %s",
                    &dim_phys_name, &tag_phys_name, name);

    if (retval != 3)
      ecs_error
        (__FILE__, __LINE__, 0,
         _("Error decoding line %ld of file\n\"%s\" :\n"
           "The description of entity <%d> physical group was expected in the form "
           "\"dimension tag names\"."),
           (long)(*line_num), ecs_file_get_name(mesh_file));

    if (dim_phys_name > 1) {
      if (tag_phys_name > _n_phys_names_max) {
        ECS_REALLOC(_phys_name_id, tag_phys_name + 1, int);
        for (int j = _n_phys_names_max+1; j < tag_phys_name + 1; j++)
          _phys_name_id[j] = -1;
        _n_phys_names_max = tag_phys_name;
      }
      _phys_name_id[tag_phys_name] = cpt_phys_names;

      /* Remove extra quotes */

      size_t s_id = 0, e_id = strlen(name);
      if (e_id > 2) {
        if (name[0] == '"' && name[e_id-1] == '"') {
          s_id += 1;
          e_id -= 1;
          name[e_id] = '\0';
        }
      }
      ECS_REALLOC(_phys_name, cpt_phys_names+1, char *);
      ECS_MALLOC(_phys_name[cpt_phys_names], e_id - s_id + 1, char);
      strcpy(_phys_name[cpt_phys_names], name+s_id);

      cpt_phys_names += 1;
    }

  }

  /* End of section line */

  ecs_file_gets(line, ECS_LOC_LNG_MAX_CHAINE_GMSH, mesh_file, line_num);

  if (strncmp(line, "$EndPhysicalNames", strlen("$EndPhysicalNames")) != 0)
    ecs_error(__FILE__, __LINE__, 0,
              _("The end of the physical names specification of file\n"
                "\"%s\" is not present at the expected place.)"),
              ecs_file_get_name(mesh_file));

  /* Return values */

  *n_phys_names = cpt_phys_names;
  *n_phys_names_max = _n_phys_names_max;
  *phys_name_id = _phys_name_id;
  *phys_name = _phys_name;
}

/*----------------------------------------------------------------------------
 * Read entities
 *----------------------------------------------------------------------------*/

static void
_read_entities(int             n_phys_names_max,
               const int       phys_name_id[],
               const char    **phys_name,
               ecs_file_t     *mesh_file,
               int            *line_num,
               int             version_fmt_gmsh,
               int             type_fmt_gmsh,
               int             fam_shift[ECS_N_ENTMAIL],
               ecs_famille_t  *vect_famille_tete[ECS_N_ENTMAIL])
{
  char   line[ECS_LOC_LNG_MAX_CHAINE_GMSH];
  int    retval;

  unsigned long  _n_elts[4] = {0, 0, 0, 0};

  int  n_vtx_b_coo = (version_fmt_gmsh < 41) ? 6 : 3;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Initialize families */

  for (int l = 0; l < ECS_N_ENTMAIL; l++) {
    vect_famille_tete[l] = NULL;
    fam_shift[l] = 0;
  }

  /* Read based on format */

  if (type_fmt_gmsh == 0) {

    ecs_file_gets(line, ECS_LOC_LNG_MAX_CHAINE_GMSH, mesh_file, line_num);

    retval = sscanf(line,"%lu %lu %lu %lu",
                    &(_n_elts[0]), &(_n_elts[1]), &(_n_elts[2]), &(_n_elts[3]));

    if (retval != 4)
      ecs_error(__FILE__, __LINE__, 0,
                _("Error decoding line %ld of file\n\"%s\" :\n"
                  "The information on entities form "
                  "\"numPoints numCurves numSurfaces numVolumes\"."),
                (long)(*line_num), ecs_file_get_name(mesh_file));

  }
  else {

    ecs_file_set_type(mesh_file, ECS_FILE_TYPE_BINARY);

    ecs_file_read(_n_elts, sizeof(unsigned long), 4, mesh_file);

  }

  /* Loop on entity type */

  int ifam_base = 0;

  for (int l = 0; l < 4; l++) {

    int  ifam_max = 0;
    ecs_int_t  ifam_ent = -1;

    const int n_coo = (l == 0) ? n_vtx_b_coo : 6;

    if (l == 2)
      ifam_ent = ECS_ENTMAIL_FAC;
    else if (l == 3)
      ifam_ent = ECS_ENTMAIL_CEL;

    if (ifam_ent > -1)
      fam_shift[ifam_ent] = ifam_base;

    unsigned int  n = _n_elts[l];

    for (unsigned long k = 0; k < n; k++) {

      int err = 0;

      /* Text variant */

      if (type_fmt_gmsh == 0) {

        ecs_file_gets(line, ECS_LOC_LNG_MAX_CHAINE_GMSH, mesh_file, line_num);

        char *ssch = strtok(line, " ");
        if (ssch == NULL) {
          err = 1;
          break;
        }

        int tag = atoi(ssch);

        if (ifam_ent > -1 && tag > ifam_max)
          ifam_max = tag;

        for (int m = 0; m < n_coo; m++) { /* coords or bounding box (unused here) */
          ssch = strtok(NULL, " ");
          if (ssch == NULL) {
            err = 1;
            break;
          }
        }

        ssch = strtok(NULL, " ");
        if (ssch == NULL) {
          err = 1;
          break;
        }

        unsigned long n_phys = atol(ssch);

        ecs_descr_t  *descr_tete = NULL;

        for (unsigned long m = 0; m < n_phys; m++) {
          ssch = strtok(NULL, " ");
          if (ssch == NULL) {
            err = 1;
            break;
          }
          int phys_tag = atoi(ssch);
          if (ifam_ent > -1) {
            ecs_descr_t *descr = ecs_descr__cree(phys_tag, NULL);
            ecs_descr_chaine__ajoute(&descr_tete, descr);
            if (phys_tag > -1 && phys_tag <= n_phys_names_max) {
              int id = phys_name_id[phys_tag];
              if (id > -1) {
                const char *name = phys_name[id];
                descr = ecs_descr__cree(ECS_DESCR_IDE_NUL, name);
                ecs_descr_chaine__ajoute(&descr_tete, descr);
              }
            }
          }
        }

        if (ifam_ent > -1) {
          ecs_famille_t *fml = ecs_famille__cree(tag + ifam_base, descr_tete);
          ecs_famille_chaine__ajoute(&vect_famille_tete[ifam_ent], fml);
        }

        if (l > 0) { /* bounding entities */

          ssch = strtok(NULL, " ");
          if (ssch == NULL) {
            err = 1;
            break;
          }

          unsigned long n_bound = atol(ssch);
          for (unsigned long m = 0; m < n_bound; m++) {
            ssch = strtok(NULL, " ");
            if (ssch == NULL) {
              err = 1;
              break;
            }
          }

        }

        if (err > 0)
          ecs_error(__FILE__, __LINE__, 0,
                    _("Error decoding line %ld of file\n\"%s\" :\n"
                      "The information on entities does not seem to match the spec."),
                    (long)(*line_num), ecs_file_get_name(mesh_file));

      }
      else { /* Binary variant */

        int tag = 0;
        ecs_file_read(&tag, sizeof(int), 1, mesh_file);

        if (ifam_ent > -1 && tag > ifam_max)
          ifam_max = tag;

        /* bounding box (unused here) */
        double bbox[6];
        ecs_file_read(bbox, sizeof(double), n_coo, mesh_file);

        unsigned long n_phys = 0;
        ecs_file_read(&n_phys, sizeof(unsigned long), 1, mesh_file);

        ecs_descr_t  *descr_tete = NULL;

        for (unsigned long m = 0; m < n_phys; m++) {
          int phys_tag = 0;
          ecs_file_read(&phys_tag, sizeof(int), 1, mesh_file);
          if (ifam_ent > -1) {
            ecs_descr_t *descr = ecs_descr__cree(phys_tag, NULL);
            ecs_descr_chaine__ajoute(&descr_tete, descr);
            if (phys_tag > -1 && phys_tag <= n_phys_names_max) {
              int id = phys_name_id[phys_tag];
              if (id > -1) {
                const char *name = phys_name[id];
                descr = ecs_descr__cree(ECS_DESCR_IDE_NUL, name);
                ecs_descr_chaine__ajoute(&descr_tete, descr);
              }
            }
          }
        }

        if (ifam_ent > -1) {
          ecs_famille_t *fml = ecs_famille__cree(tag + ifam_base, descr_tete);
          ecs_famille_chaine__ajoute(&vect_famille_tete[ifam_ent], fml);
        }

        if (l > 0) { /* bounding entities */

          unsigned long n_bound = 0;
          ecs_file_read(&n_bound, sizeof(unsigned long), 1, mesh_file);

          for (unsigned long m = 0; m < n_bound; m++) {
            int bound_tag;
            ecs_file_read(&bound_tag, sizeof(int), 1, mesh_file);
          }

        }

      } /* End of test on format */

    } /* End of loop on elements */

    if (ifam_ent > -1) {
      int j = 0;
      ifam_max += ifam_base;
      while (ifam_max > 0) {
        j++;
        ifam_max /= 10;
      }
      ifam_max = 1;
      while (j > 0) {
        ifam_max *= 10;
        j--;
      }
      ifam_base = ifam_max;
    }

  } /* end of loop on entity types */

  /* End of section line */

  if (type_fmt_gmsh != 0) {
    ecs_file_set_type(mesh_file, ECS_FILE_TYPE_TEXT);
    ecs_file_gets(line, 2, mesh_file, line_num);
  }

  ecs_file_gets(line, ECS_LOC_LNG_MAX_CHAINE_GMSH, mesh_file, line_num);

  if (strncmp(line, "$EndEntities", strlen("$EndEntities")) != 0)
    ecs_error(__FILE__, __LINE__, 0,
              _("The end of the entities specification of file\n"
                "\"%s\" is not present at the expected place.)"),
              ecs_file_get_name(mesh_file));
}

/*----------------------------------------------------------------------------
 *  Lecture des coordonnées des noeuds
 *----------------------------------------------------------------------------*/

static void
ecs_loc_pre_gmsh__lit_nodes_v4(ecs_maillage_t   *maillage,
                               int               version_fmt_gmsh,
                               ecs_file_t       *fic_maillage,
                               int               type_fmt_gmsh,
                               int              *num_ligne,
                               ecs_int_t       **som_val_label)
{
  char   chaine[ECS_LOC_LNG_MAX_CHAINE_GMSH];
  int    retour;

  /* Variables Gmsh lues */

  double         coord[3];
  unsigned long  n_ent_blocks, nbr_nod;

  ecs_int_t      icoo;

  /* Stockage avant transfert */

  ecs_int_t    ind_nod = 0;

  ecs_coord_t * som_val_coord;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Initialisations */
  /*=================*/

  /* Boucle de lecture des noeuds */
  /*==============================*/

  if (type_fmt_gmsh == 0) {

    ecs_file_gets(chaine, ECS_LOC_LNG_MAX_CHAINE_GMSH, fic_maillage, num_ligne);

    retour = sscanf(chaine,"%lu %lu", &n_ent_blocks, &nbr_nod);

    ECS_MALLOC(som_val_coord, nbr_nod * 3, ecs_coord_t);
    ECS_MALLOC(*som_val_label, nbr_nod, ecs_int_t);

    if (retour != 2)
      ecs_error(__FILE__, __LINE__, 0,
                _("Error decoding line %ld of file\n\"%s\" :\n"
                  "The information on node counts was expected in the form "
                  "\"numEntityBlocks numNodes\"."),
                (long)(*num_ligne), ecs_file_get_name(fic_maillage));

    for (unsigned long eb = 0; eb < n_ent_blocks; eb++) {

      int tag_ent, dim_ent, is_param;
      unsigned long n_ent_nodes;

      ecs_file_gets(chaine,
                    ECS_LOC_LNG_MAX_CHAINE_GMSH,
                    fic_maillage,
                    num_ligne);

      retour = sscanf(chaine,"%d %d %d %lu",
                      &tag_ent, &dim_ent, &is_param, &n_ent_nodes);

      if (retour != 4)
        ecs_error
          (__FILE__, __LINE__, 0,
           _("Error decoding line %ld of file\n\"%s\" :\n"
             "The description of entity <%d> points was expected in the form "
             "\"tagEntity dimEntity parametric numNodes\"."),
           (long)(*num_ligne), ecs_file_get_name(fic_maillage),
           eb);

      if (version_fmt_gmsh < 41) {

        for (unsigned long bn = 0; bn < n_ent_nodes; bn++, ind_nod++) {

          ecs_file_gets(chaine,
                        ECS_LOC_LNG_MAX_CHAINE_GMSH,
                        fic_maillage,
                        num_ligne);

          int label;
          retour = sscanf(chaine,"%d %lg %lg %lg",
                          &label, coord, coord+1, coord+2);

          if (retour != 4)
            ecs_error(__FILE__, __LINE__, 0,
                      _("Error decoding line %ld of file\n\"%s\" :\n"
                        "The description of point <%ld> was expected in the form "
                        "\"label x y z\"."),
                      (long)(*num_ligne), ecs_file_get_name(fic_maillage),
                      (long)(ind_nod+1));

          /* Étiquette du noeud lu */

          (*som_val_label)[ind_nod] = (ecs_int_t)label;

          /* Coordonnées du noeud lu */

          for (icoo = 0; icoo < 3; icoo++)
            som_val_coord[ind_nod * 3 + icoo] = (ecs_coord_t)(coord[icoo]);

        }

      }
      else {

        for (unsigned long bn = 0; bn < n_ent_nodes; bn++, ind_nod++) {

          ecs_file_gets(chaine,
                        ECS_LOC_LNG_MAX_CHAINE_GMSH,
                        fic_maillage,
                        num_ligne);

          unsigned long label;
          retour = sscanf(chaine, "%lu", &label);

          if (retour != 1)
            ecs_error(__FILE__, __LINE__, 0,
                      _("Error decoding line %ld of file\n\"%s\" :\n"
                        "The point <%ld> label was expected in the form "
                        "\"label\"."),
                      (long)(*num_ligne), ecs_file_get_name(fic_maillage),
                      (long)(ind_nod+1));

          /* Étiquette du noeud lu */

          (*som_val_label)[ind_nod] = (ecs_int_t)label;

        }

        ind_nod -= n_ent_nodes;

        for (unsigned long bn = 0; bn < n_ent_nodes; bn++, ind_nod++) {

          ecs_file_gets(chaine,
                        ECS_LOC_LNG_MAX_CHAINE_GMSH,
                        fic_maillage,
                        num_ligne);

          retour = sscanf(chaine, "%lg %lg %lg", coord, coord+1, coord+2);

          if (retour != 3)
            ecs_error(__FILE__, __LINE__, 0,
                      _("Error decoding line %ld of file\n\"%s\" :\n"
                        "The point <%ld> coordinates were expected in the form "
                        "\"x y z\"."),
                      (long)(*num_ligne), ecs_file_get_name(fic_maillage),
                      (long)(ind_nod+1));

          /* Coordonnées du noeud lu */

          for (icoo = 0; icoo < 3; icoo++)
            som_val_coord[ind_nod * 3 + icoo] = (ecs_coord_t)(coord[icoo]);

        }

      }

    }

  }
  else {

    ecs_file_set_type(fic_maillage, ECS_FILE_TYPE_BINARY);

    if (version_fmt_gmsh < 41) {
      unsigned long buf[2];
      ecs_file_read(buf, sizeof(unsigned long), 2, fic_maillage);

      n_ent_blocks = buf[0];
      nbr_nod = buf[1];
    }
    else {
      size_t buf[4];
      ecs_file_read(buf, sizeof(size_t), 4, fic_maillage);

      n_ent_blocks = buf[0];
      nbr_nod = buf[1];
    }

    ECS_MALLOC(som_val_coord, nbr_nod * 3, ecs_coord_t);
    ECS_MALLOC(*som_val_label, nbr_nod, ecs_int_t);

    for (unsigned long eb = 0; eb < n_ent_blocks; eb++) {

      int buf_i[3];

      ecs_file_read(buf_i, sizeof(int), 3, fic_maillage);

      if (version_fmt_gmsh < 41) {

        unsigned long n_ent_nodes;
        ecs_file_read(&n_ent_nodes, sizeof(unsigned long), 1, fic_maillage);

        for (unsigned long bn = 0; bn < n_ent_nodes; bn++, ind_nod++) {

          int label;
          double xyz[3];

          ecs_file_read(&label, sizeof(int), 1, fic_maillage);
          ecs_file_read(&xyz, sizeof(double), 3, fic_maillage);

          (*som_val_label)[ind_nod] = label;
          for (icoo = 0; icoo < 3; icoo++)
            som_val_coord[ind_nod * 3 + icoo] = xyz[icoo];

        }

      }
      else {

        size_t n_ent_nodes;
        ecs_file_read(&n_ent_nodes, sizeof(size_t), 1, fic_maillage);

        for (unsigned long bn = 0; bn < n_ent_nodes; bn++, ind_nod++) {
          size_t label;
          ecs_file_read(&label, sizeof(size_t), 1, fic_maillage);
          (*som_val_label)[ind_nod] = label;
        }

        ind_nod -= n_ent_nodes;

        for (unsigned long bn = 0; bn < n_ent_nodes; bn++, ind_nod++) {
          double xyz[3];
          ecs_file_read(&xyz, sizeof(double), 3, fic_maillage);
          for (icoo = 0; icoo < 3; icoo++)
            som_val_coord[ind_nod * 3 + icoo] = xyz[icoo];
        }

      }

    }

    ecs_file_set_type(fic_maillage, ECS_FILE_TYPE_TEXT);

    ecs_file_gets(chaine, 2, fic_maillage, num_ligne);

  }

  /* Ligne de fin de rubrique */

  ecs_file_gets(chaine, ECS_LOC_LNG_MAX_CHAINE_GMSH, fic_maillage, num_ligne);

  if (   (strncmp(chaine, "$EndNodes", strlen("$EndNodes")) != 0)
      && (strncmp(chaine, "$ENDNOD", strlen("$ENDNOD")) != 0))
    ecs_error(__FILE__, __LINE__, 0,
              _("The end of the points list specification of file\n"
                "\"%s\" is not present at the expected place.)"),
              ecs_file_get_name(fic_maillage));

  /* Transfert des valeurs lues dans la structure d'entité de maillage */
  /*===================================================================*/

  ecs_maillage_pre__cree_som(maillage,
                             nbr_nod,
                             som_val_coord);
}

/*----------------------------------------------------------------------------
 *  Lecture de la table de connectivité
 *----------------------------------------------------------------------------*/

static void
ecs_loc_pre_gmsh__lit_elements(ecs_maillage_t  *maillage,
                               ecs_file_t      *fic_maillage,
                               int              version_fmt_gmsh,
                               int              type_fmt_gmsh,
                               int             *num_ligne,
                               ecs_int_t        nbr_som,
                               ecs_int_t      **som_val_label)
{
  char        chaine[ECS_LOC_LNG_MAX_CHAINE_GMSH];
  char       *ssch;

  /* Variables Gmsh lues */

  long        label;
  ecs_int_t   type_ecs;
  ecs_int_t   coul_elt;
  ecs_int_t   ind_elt;
  ecs_int_t   nbr_elt;
  ecs_int_t   nbr_elt_lus;
  ecs_int_t   nbr_elt_tot;
  ecs_int_t   ind_nod_elt;
  ecs_int_t   ind_som_elt;
  ecs_int_t   nbr_som_elt;
  ecs_int_t   ind_tag_elt;
  ecs_int_t   nbr_tag_elt;

  gmsh_int_t  header[3];
  gmsh_int_t  data[100];

  int         type_gmsh;
  ecs_int_t   nbr_nod_elt_gmsh;
  ecs_int_t   num_nod_elt_gmsh[ECS_GMSH_NBR_MAX_SOM];
  ecs_int_t   ent_num;

  ecs_int_t   icoul;
  ecs_int_t   ient;

  ecs_int_t    cpt_coul_ent[ECS_N_ENTMAIL]; /* Compteur de couleurs         */
  ecs_int_t   *val_coul_ent[ECS_N_ENTMAIL]; /* Tableau valeurs des couleurs */
  ecs_size_t  *cpt_elt_coul_ent[ECS_N_ENTMAIL];

  /* Stockage avant transfert */
  /*--------------------------*/

  size_t       cpt_elt_ent        [ECS_N_ENTMAIL]; /* Nombre d'elems/entite */

  ecs_size_t  *elt_pos_som_ent    [ECS_N_ENTMAIL]; /* Positions numeros som */
  ecs_int_t   *elt_val_som_ent    [ECS_N_ENTMAIL]; /* Numeros des sommets   */
  ecs_int_t   *elt_val_color_ent  [ECS_N_ENTMAIL]; /* Couleurs              */

  ecs_int_t cpt_are   = 0;
  ecs_int_t cpt_point = 0;

  bool       ligne_decodee = true;

  int _version_fmt = version_fmt_gmsh / 10;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Initialisations */
  /*=================*/

  /* Décodage du nombre d'éléments (toutes dimensions confondues */

  ecs_file_gets(chaine, ECS_LOC_LNG_MAX_CHAINE_GMSH, fic_maillage, num_ligne);

  nbr_elt = (ecs_int_t)(atol(chaine));

  /* Allocations des tableaux locaux */
  /*=================================*/

  /* Attention au decalage de `1' !!!         */
  /* On n'alloue pas les tableaux locaux pour */
  /* `ECS_ENTMAIL_DEB = ECS_ENTMAIL_SOM'      */

  for (ient = ECS_ENTMAIL_FAC; ient < ECS_N_ENTMAIL; ient++) {

    cpt_elt_ent        [ient] = 0;

    cpt_coul_ent       [ient] = 0   ;
    cpt_elt_coul_ent   [ient] = NULL;
    val_coul_ent       [ient] = NULL;

  }

  /*
    On surdimensionne les tableaux de lecture des éléments; on les
    redimensionnera à la fin.
  */

  for (ient = ECS_ENTMAIL_FAC; ient < ECS_N_ENTMAIL; ient++) {

    ECS_MALLOC(elt_val_color_ent[ient],   nbr_elt, ecs_int_t);
    ECS_MALLOC(elt_pos_som_ent[ient],     nbr_elt + 1, ecs_size_t);

    elt_pos_som_ent[ient][0] = 1;
  }

  /*
    Plus "gros" élément linéaire disponible :
    - quadrangle (4 sommets) pour les faces
    - hexaèdre (8 sommets) pour les cellules.
  */

  ECS_MALLOC(elt_val_som_ent[ECS_ENTMAIL_FAC], nbr_elt * 4, ecs_int_t);
  ECS_MALLOC(elt_val_som_ent[ECS_ENTMAIL_CEL], nbr_elt * 8, ecs_int_t);


  /* Boucle de lecture des éléments */
  /*================================*/

  if (type_fmt_gmsh == 0) {

    for (ind_elt = 0; ind_elt < nbr_elt; ind_elt++) {

      if (ligne_decodee == false)
        break;

      ligne_decodee = false;

      ecs_file_gets(chaine, ECS_LOC_LNG_MAX_CHAINE_GMSH,
                    fic_maillage, num_ligne);

      /* Au format Gmsh 1.0, pour chaque élément, on a les entiers suivants :
         label type num_ent_phys num_ent_elem nb_sommets <liste_sommets> */

      /* Au format Gmsh 2.0, pour chaque élément, on a les entiers suivants :
         label type nb_tags <tag> <liste_sommets> */

      /* Lecture du label de l'élément courant */

      ssch = strtok(chaine, " ");
      if (ssch == NULL)
        break;

      label = atol(ssch);

      /* Lecture du type de l'élément courant */

      ssch   = strtok(NULL, " ");
      if (ssch == NULL)
        break;

      type_gmsh = atoi(ssch);

      if (type_gmsh < (int)GMSH_SEG2 || type_gmsh > (int) GMSH_POINT1)
        ecs_error(__FILE__, __LINE__, 0,
                  _("Error reading a Gmsh mesh file:\n"
                    "at line %ld of file \"%s\".\n"
                    "Type identifier <%d> for element <%ld> is not recognized."),
                  (long)(*num_ligne), ecs_file_get_name(fic_maillage),
                  (int)type_gmsh, (long)label);

      type_ecs         = ecs_gmsh_elt_liste_c[type_gmsh - 1].ecs_typ;
      nbr_nod_elt_gmsh = ecs_gmsh_elt_liste_c[type_gmsh - 1].nbr_som;
      nbr_som_elt      = ecs_fic_elt_typ_liste_c[type_ecs].nbr_som;

      if (type_gmsh == GMSH_POINT1) {
        cpt_point += 1;
        ligne_decodee = true;
        continue;
      }
      else if (type_gmsh == GMSH_SEG2 || type_gmsh == GMSH_SEG3) {
        cpt_are += 1;
        ligne_decodee = true;
        continue;
      }

      /* Lecture des "tags" de l'élément courant */

      coul_elt = 0;

      if (_version_fmt == 2) {

        ssch   = strtok(NULL, " ");
        if (ssch == NULL)
          break;

        nbr_tag_elt = (ecs_int_t)(atol(ssch));

        for (ind_tag_elt = 0; ind_tag_elt < nbr_tag_elt; ind_tag_elt++) {

          ssch   = strtok(NULL, " ");
          if (ssch == NULL)
            break;

          /* Par défaut, Gmsh écrit des fichiers avec 2 "tags",
             le premier correspondant à un numéro d'entité physique,
             le second à une entité géométrique élémentaire;
             on associe le premier à une couleur */

          if (ind_tag_elt == 0)
            coul_elt = (ecs_int_t)(atol(ssch));

        }

        if (ssch == NULL && ind_tag_elt < nbr_tag_elt)
          break;

      }
      else { /* Ancienne version 1 du format Gmsh */

        /* On interprête reg-phys (numéro d'entité physique) comme
           étant une couleur */

        ssch   = strtok(NULL, " ");
        if (ssch == NULL)
          break;

        coul_elt = (ecs_int_t)(atol(ssch));

        /* On saute la valeur reg-elem (numéro d'entité élémentaire) */

        ssch   = strtok(NULL, " ");
        if (ssch == NULL)
          break;

        /* On saute le nombre de noeuds, redondant */

        ssch   = strtok(NULL, " ");
        if (ssch == NULL)
          break;

      }

      /* Lecture des numéros des sommets de l'élément courant */

      for (ind_nod_elt = 0; ind_nod_elt < nbr_nod_elt_gmsh; ind_nod_elt++) {

        ssch   = strtok(NULL, " ");
        if (ssch == NULL)
          break;

        num_nod_elt_gmsh[ind_nod_elt] = (ecs_int_t)(atol(ssch));
      }

      if (ssch == NULL && ind_nod_elt < nbr_nod_elt_gmsh)
        break;

      ligne_decodee = true;

      /* Stockage des valeurs avant transfert dans la structure `maillage' */
      /*===================================================================*/

      /* Identification de l'entité concernée */

      ent_num = ecs_maillage_pre__ret_typ_geo(type_ecs);

      /* Position des numéros de sommets du prochain élément */

      elt_pos_som_ent[ent_num][cpt_elt_ent[ent_num] + 1] =
        elt_pos_som_ent[ent_num][cpt_elt_ent[ent_num]] + nbr_som_elt;

      /* Connectivité de l'élément par ses numéros de sommets */

      for (ind_som_elt = 0; ind_som_elt < nbr_som_elt; ind_som_elt++) {

        elt_val_som_ent
          [ent_num]
          [elt_pos_som_ent[ent_num][cpt_elt_ent[ent_num]] - 1 + ind_som_elt]
          = num_nod_elt_gmsh
          [ecs_gmsh_elt_liste_c[type_gmsh - 1].num_som[ind_som_elt] - 1];

      }

      /* Couleur (tag) de l'élément lu */

      icoul = 0;
      while (icoul < cpt_coul_ent[ent_num]           &&
             val_coul_ent[ent_num][icoul] != coul_elt)
        icoul++;

      if (icoul == cpt_coul_ent[ent_num]) {

        /* La valeur de la couleur n'a pas encore été stockée */

        ECS_REALLOC(val_coul_ent[ent_num]    , cpt_coul_ent[ent_num] + 1,
                    ecs_int_t);
        ECS_REALLOC(cpt_elt_coul_ent[ent_num], cpt_coul_ent[ent_num] + 1,
                    ecs_size_t);
        cpt_elt_coul_ent[ent_num][icoul] = 0;
        val_coul_ent[ent_num][icoul] = coul_elt;
        cpt_coul_ent[ent_num]++;

      }

      cpt_elt_coul_ent[ent_num][icoul]++;
      elt_val_color_ent[ent_num][cpt_elt_ent[ent_num]] = icoul + 1;

      /* Incrémentation du nombre d'éléments lus */

      cpt_elt_ent[ent_num]++;

    }

    /* Fin de la boucle de lecture sur les éléments
       (et sur les lignes du fichier) */

    if (ligne_decodee == false)
      ecs_error(__FILE__, __LINE__, 0,
                _("Error decoding line %ld "
                  "(corresponding to element <%ld>) of file\n"
                  "\"%s\"."),
                (long)(*num_ligne), (long)(ind_elt+1),
                ecs_file_get_name(fic_maillage));

  }
  else {

    ecs_file_set_type(fic_maillage, ECS_FILE_TYPE_BINARY);

    nbr_elt_lus = 0;
    nbr_elt_tot = nbr_elt;

    while (nbr_elt_lus != nbr_elt_tot) {

      ecs_file_read(&header, sizeof(gmsh_int_t), 3, fic_maillage);

      type_gmsh   = header[0];
      nbr_elt     = header[1];
      nbr_tag_elt = header[2];

      if (type_gmsh < (int)GMSH_SEG2 || type_gmsh > (int) GMSH_POINT1)
        ecs_error(__FILE__, __LINE__, 0,
                  _("Error reading a Gmsh mesh file:\n"
                    "Type identifier <%d> is not recognized."),
                  ecs_file_get_name(fic_maillage), (int)type_gmsh);

      type_ecs         = ecs_gmsh_elt_liste_c[type_gmsh - 1].ecs_typ;
      nbr_nod_elt_gmsh = ecs_gmsh_elt_liste_c[type_gmsh - 1].nbr_som;
      nbr_som_elt      = ecs_fic_elt_typ_liste_c[type_ecs].nbr_som;

      for (ind_elt = 0; ind_elt < nbr_elt; ind_elt++) {

        ecs_file_read(&data, sizeof(gmsh_int_t),
                      1 + nbr_tag_elt + nbr_nod_elt_gmsh, fic_maillage);

        if (type_gmsh == GMSH_POINT1) {
          cpt_point += 1;
          continue;
        }
        else if (type_gmsh == GMSH_SEG2 || type_gmsh == GMSH_SEG3) {
          cpt_are += 1;
          continue;
        }

        /* Par défaut, Gmsh écrit des fichiers avec 2 "tags",
           le premier correspondant à un numéro d'entité physique,
           le second à une entité géométrique élémentaire;
           on associe le premier à une couleur */

        label    = data[0];
        coul_elt = data[1];

        /* Lecture des numéros des sommets de l'élément courant */

        for (ind_nod_elt = 0; ind_nod_elt < nbr_nod_elt_gmsh; ind_nod_elt++)
          num_nod_elt_gmsh[ind_nod_elt] = data[1 + nbr_tag_elt + ind_nod_elt];

        /* Stockage des valeurs avant transfert dans la structure `maillage' */
        /*===================================================================*/

        /* Identification de l'entité concernée */

        ent_num = ecs_maillage_pre__ret_typ_geo(type_ecs);

        /* Position des numéros de sommets du prochain élément */

        elt_pos_som_ent[ent_num][cpt_elt_ent[ent_num] + 1] =
          elt_pos_som_ent[ent_num][cpt_elt_ent[ent_num]] + nbr_som_elt;

        /* Connectivité de l'élément par ses numéros de sommets */

        for (ind_som_elt = 0; ind_som_elt < nbr_som_elt; ind_som_elt++) {

          elt_val_som_ent
            [ent_num]
            [elt_pos_som_ent[ent_num][cpt_elt_ent[ent_num]] - 1 + ind_som_elt]
            = num_nod_elt_gmsh
            [ecs_gmsh_elt_liste_c[type_gmsh - 1].num_som[ind_som_elt] - 1];

        }

        /* Couleur (tag) de l'élément lu */

        icoul = 0;
        while (icoul < cpt_coul_ent[ent_num]           &&
               val_coul_ent[ent_num][icoul] != coul_elt)
          icoul++;

        if (icoul == cpt_coul_ent[ent_num]) {

          /* La valeur de la couleur n'a pas encore été stockée */

          ECS_REALLOC(val_coul_ent[ent_num]    , cpt_coul_ent[ent_num] + 1,
                      ecs_int_t);
          ECS_REALLOC(cpt_elt_coul_ent[ent_num], cpt_coul_ent[ent_num] + 1,
                      ecs_size_t);
          cpt_elt_coul_ent[ent_num][icoul] = 0;
          val_coul_ent[ent_num][icoul] = coul_elt;
          cpt_coul_ent[ent_num]++;

        }

        cpt_elt_coul_ent[ent_num][icoul]++;
        elt_val_color_ent[ent_num][cpt_elt_ent[ent_num]] = icoul + 1;

        /* Incrémentation du nombre d'éléments lus */

        cpt_elt_ent[ent_num]++;

      }

      nbr_elt_lus += nbr_elt;
    }

    ecs_file_set_type(fic_maillage, ECS_FILE_TYPE_TEXT);

    ecs_file_gets(chaine, 2, fic_maillage, num_ligne);

  }

  /* Indications */

  if (cpt_point > -1 || cpt_are > 0)
    printf("\n");
  if (cpt_point > -1)
    printf(_("  %ld elements of type \"point\" ignored\n"),
           (long)cpt_point);
  if (cpt_are > -1)
    printf(_("  %ld elements of type \"edge\" ignored\n"),
           (long)cpt_are);
  if (cpt_point > -1 || cpt_are > 0)
    printf("\n");

  /* Ligne de fin de rubrique */

  ecs_file_gets(chaine, ECS_LOC_LNG_MAX_CHAINE_GMSH, fic_maillage, num_ligne);

  if (   (strncmp(chaine, "$EndElements", strlen("$EndElements")) != 0)
      && (strncmp(chaine, "$ENDELM", strlen("$ENDELM")) != 0))
    ecs_error(__FILE__, __LINE__, 0,
              _("The end of the elements list specification of file\n"
                "\"%s\" is not present at the expected place.)"),
              ecs_file_get_name(fic_maillage));

  /* Réallocations des tableaux locaux */
  /*===================================*/

  /* Reallocations des tableaux locaux et mise à jour des références */
  /*=================================================================*/

  for (ient = ECS_ENTMAIL_FAC; ient < ECS_N_ENTMAIL; ient++) {

    if (cpt_elt_ent[ient] != 0) {
      ECS_REALLOC(elt_pos_som_ent[ient]    ,
                  cpt_elt_ent[ient] + 1    , ecs_size_t);
      ECS_REALLOC(elt_val_som_ent[ient]    ,
                  elt_pos_som_ent[ient][cpt_elt_ent[ient]] - 1, ecs_int_t);
      ECS_REALLOC(elt_val_color_ent[ient]  ,
                  cpt_elt_ent[ient]        , ecs_int_t);
    }

    ecs_maillage_pre__label_en_indice
      (nbr_som,
       elt_pos_som_ent[ient][cpt_elt_ent[ient]] - 1,
       *som_val_label,
       elt_val_som_ent[ient]);
  }

  ECS_FREE(*som_val_label);

  /* Transfert des valeurs lues dans les structures d'entités de maillage */
  /*======================================================================*/

  ecs_maillage_pre__cree_elt(maillage,
                             cpt_elt_ent,
                             elt_pos_som_ent,
                             elt_val_som_ent,
                             NULL,
                             elt_val_color_ent,
                             cpt_coul_ent,
                             val_coul_ent,
                             cpt_elt_coul_ent);
}

/*----------------------------------------------------------------------------
 *  Lecture de la table de connectivité
 *----------------------------------------------------------------------------*/

static void
ecs_loc_pre_gmsh__lit_elements_v4(ecs_maillage_t  *maillage,
                                  ecs_file_t      *fic_maillage,
                                  int              version_fmt_gmsh,
                                  int              type_fmt_gmsh,
                                  int             *num_ligne,
                                  ecs_int_t        nbr_som,
                                  int              fam_shift[ECS_N_ENTMAIL],
                                  ecs_int_t      **som_val_label)
{
  char        chaine[ECS_LOC_LNG_MAX_CHAINE_GMSH];
  char       *ssch;

  /* Variables Gmsh lues */

  unsigned long  nbr_elt;
  unsigned long  n_ent_blocks;

  long        label;
  ecs_int_t   type_ecs;
  ecs_int_t   nbr_elt_lus;
  ecs_int_t   ind_nod_elt;
  ecs_int_t   ind_som_elt;
  ecs_int_t   nbr_som_elt;

  int         type_gmsh;
  ecs_int_t   nbr_nod_elt_gmsh;
  ecs_int_t   num_nod_elt_gmsh[ECS_GMSH_NBR_MAX_SOM];
  ecs_int_t   ent_num;

  ecs_int_t   ient;

  /* Stockage avant transfert */
  /*--------------------------*/

  size_t       cpt_elt_ent        [ECS_N_ENTMAIL]; /* Nombre d'elems/entite */

  ecs_size_t  *elt_pos_som_ent    [ECS_N_ENTMAIL]; /* Positions numeros som */
  ecs_int_t   *elt_val_som_ent    [ECS_N_ENTMAIL]; /* Numeros des sommets   */
  int         *elt_val_fam_ent    [ECS_N_ENTMAIL]; /* Groupes (familles)    */

  ecs_int_t  ind_elt = 0;
  ecs_int_t  cpt_are   = 0;
  ecs_int_t  cpt_point = 0;

  int        retour;
  bool       ligne_decodee = true;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Allocations des tableaux locaux */
  /*=================================*/

  /* Attention au decalage de `1' !!!         */
  /* On n'alloue pas les tableaux locaux pour */
  /* `ECS_ENTMAIL_DEB = ECS_ENTMAIL_SOM'      */

  for (ient = ECS_ENTMAIL_FAC; ient < ECS_N_ENTMAIL; ient++) {
    cpt_elt_ent[ient] = 0;
  }

  /* Décodage du nombre d'éléments (toutes dimensions confondues */

  if (type_fmt_gmsh == 0) {

    ecs_file_gets(chaine, ECS_LOC_LNG_MAX_CHAINE_GMSH, fic_maillage, num_ligne);

    retour = sscanf(chaine,"%lu %lu", &n_ent_blocks, &nbr_elt);

    if (retour != 2)
      ecs_error(__FILE__, __LINE__, 0,
                _("Error decoding line %ld of file\n\"%s\" :\n"
                  "The information on element counts was expected in the form "
                  "\"numEntityBlocks numElements\"."),
                (long)(*num_ligne), ecs_file_get_name(fic_maillage));

  }
  else {

    ecs_file_set_type(fic_maillage, ECS_FILE_TYPE_BINARY);

    if (version_fmt_gmsh < 41) {
      unsigned long buf[2];
      ecs_file_read(buf, sizeof(unsigned long), 2, fic_maillage);
      n_ent_blocks = buf[0];
      nbr_elt = buf[1];
    }
    else {
      size_t buf[4];
      ecs_file_read(buf, sizeof(size_t), 4, fic_maillage);
      n_ent_blocks = buf[0];
      nbr_elt = buf[1];
    }

  }

  /*
    On surdimensionne les tableaux de lecture des éléments; on les
    redimensionnera à la fin.
  */

  for (ient = ECS_ENTMAIL_FAC; ient < ECS_N_ENTMAIL; ient++) {

    ECS_MALLOC(elt_val_fam_ent[ient], nbr_elt, int);
    ECS_MALLOC(elt_pos_som_ent[ient], nbr_elt + 1, ecs_size_t);

    elt_pos_som_ent[ient][0] = 1;
  }

  /*
    Plus "gros" élément linéaire disponible :
    - quadrangle (4 sommets) pour les faces
    - hexaèdre (8 sommets) pour les cellules.
  */

  ECS_MALLOC(elt_val_som_ent[ECS_ENTMAIL_FAC], nbr_elt * 4, ecs_int_t);
  ECS_MALLOC(elt_val_som_ent[ECS_ENTMAIL_CEL], nbr_elt * 8, ecs_int_t);


  /* Boucle de lecture des éléments */
  /*================================*/

  if (type_fmt_gmsh == 0) {

    for (unsigned long eb = 0; eb < n_ent_blocks; eb++) {

      int tag_ent, dim_ent, elt_type;
      unsigned long n_ent_elt;

      ecs_file_gets(chaine,
                    ECS_LOC_LNG_MAX_CHAINE_GMSH,
                    fic_maillage,
                    num_ligne);

      if (version_fmt_gmsh < 41)
        retour = sscanf(chaine,"%d %d %d %lu",
                        &tag_ent, &dim_ent, &elt_type, &n_ent_elt);
      else
        retour = sscanf(chaine,"%d %d %d %lu",
                        &dim_ent, &tag_ent, &elt_type, &n_ent_elt);

      if (retour != 4)
        ecs_error
          (__FILE__, __LINE__, 0,
           _("Error decoding line %ld of file\n\"%s\" :\n"
             "The description of entity <%d> elements was expected in the form "
             "\"tagEntity dimEntity typeEle numElements\"."),
           (long)(*num_ligne), ecs_file_get_name(fic_maillage),
           eb);

      type_gmsh = elt_type;

      if (type_gmsh < (int)GMSH_SEG2 || type_gmsh > (int) GMSH_POINT1)
        ecs_error(__FILE__, __LINE__, 0,
                  _("Error reading a Gmsh mesh file:\n"
                    "at line %ld of file \"%s\".\n"
                    "Type identifier <%d> for element <%ld> is not recognized."),
                  (long)(*num_ligne), ecs_file_get_name(fic_maillage),
                  (int)type_gmsh, (long)label);

      type_ecs         = ecs_gmsh_elt_liste_c[type_gmsh - 1].ecs_typ;
      nbr_nod_elt_gmsh = ecs_gmsh_elt_liste_c[type_gmsh - 1].nbr_som;
      nbr_som_elt      = ecs_fic_elt_typ_liste_c[type_ecs].nbr_som;

      for (unsigned long be = 0; be < n_ent_elt; be++, ind_elt++) {

        if (ligne_decodee == false)
          break;

        ligne_decodee = false;

        ecs_file_gets(chaine, ECS_LOC_LNG_MAX_CHAINE_GMSH,
                      fic_maillage, num_ligne);

        /* Au format Gmsh 4.0, pour chaque élément, on a les entiers suivants :
           <tag> <liste_sommets> */

        /* Lecture du label de l'élément courant */

        ssch = strtok(chaine, " ");
        if (ssch == NULL)
          break;

        if (type_gmsh == GMSH_POINT1) {
          cpt_point += 1;
          ligne_decodee = true;
          continue;
        }
        else if (type_gmsh == GMSH_SEG2 || type_gmsh == GMSH_SEG3) {
          cpt_are += 1;
          ligne_decodee = true;
          continue;
        }

        label = atol(ssch);

        /* Lecture des numéros des sommets de l'élément courant */

        for (ind_nod_elt = 0; ind_nod_elt < nbr_nod_elt_gmsh; ind_nod_elt++) {

          ssch   = strtok(NULL, " ");
          if (ssch == NULL)
            break;

          num_nod_elt_gmsh[ind_nod_elt] = (ecs_int_t)(atol(ssch));
        }

        if (ssch == NULL && ind_nod_elt < nbr_nod_elt_gmsh)
          break;

        ligne_decodee = true;

        /* Stockage des valeurs avant transfert dans la structure `maillage' */
        /*===================================================================*/

        /* Identification de l'entité concernée */

        ent_num = ecs_maillage_pre__ret_typ_geo(type_ecs);

        /* Position des numéros de sommets du prochain élément */

        elt_pos_som_ent[ent_num][cpt_elt_ent[ent_num] + 1] =
          elt_pos_som_ent[ent_num][cpt_elt_ent[ent_num]] + nbr_som_elt;

        /* Connectivité de l'élément par ses numéros de sommets */

        for (ind_som_elt = 0; ind_som_elt < nbr_som_elt; ind_som_elt++) {

          elt_val_som_ent
            [ent_num]
            [elt_pos_som_ent[ent_num][cpt_elt_ent[ent_num]] - 1 + ind_som_elt]
            = num_nod_elt_gmsh
            [ecs_gmsh_elt_liste_c[type_gmsh - 1].num_som[ind_som_elt] - 1];

        }

        /* Couleur (tag) de l'élément lu */

        elt_val_fam_ent[ent_num][cpt_elt_ent[ent_num]]
          = tag_ent + fam_shift[ent_num];

        /* Incrémentation du nombre d'éléments lus */

        cpt_elt_ent[ent_num]++;

        /* Fin de la boucle de lecture sur les éléments
           (et sur les lignes du fichier) */

        if (ligne_decodee == false)
          ecs_error(__FILE__, __LINE__, 0,
                    _("Error decoding line %ld "
                      "(corresponding to element <%ld>) of file\n"
                      "\"%s\"."),
                    (long)(*num_ligne), (long)(ind_elt+1),
                    ecs_file_get_name(fic_maillage));

      }

    }

  }
  else {

    nbr_elt_lus = 0;

    for (unsigned long eb = 0; eb < n_ent_blocks; eb++) {

      int tag_ent = -1;
      unsigned long n_ent_elt = 0;

      int buf_i[3];

      ecs_file_read(buf_i, sizeof(int), 3, fic_maillage);

      if (version_fmt_gmsh < 41) {
        tag_ent = buf_i[0];
        /* dim_ent = buf_i[1] */
      }
      else {
        /* dim_ent = buf_i[0] */
        tag_ent = buf_i[1];
      }
      type_gmsh  = buf_i[2];

      if (version_fmt_gmsh < 41) {
        ecs_file_read(&n_ent_elt, sizeof(unsigned long), 1, fic_maillage);
      }
      else {
        size_t _n_ent_elt;
        ecs_file_read(&_n_ent_elt, sizeof(size_t), 1, fic_maillage);
        n_ent_elt = _n_ent_elt;
      }

      if (type_gmsh < (int)GMSH_SEG2 || type_gmsh > (int) GMSH_POINT1)
        ecs_error(__FILE__, __LINE__, 0,
                  _("Error reading a Gmsh mesh file:\n"
                    "Type identifier <%d> is not recognized."),
                  ecs_file_get_name(fic_maillage), (int)type_gmsh);

      type_ecs         = ecs_gmsh_elt_liste_c[type_gmsh - 1].ecs_typ;
      nbr_nod_elt_gmsh = ecs_gmsh_elt_liste_c[type_gmsh - 1].nbr_som;
      nbr_som_elt      = ecs_fic_elt_typ_liste_c[type_ecs].nbr_som;

      for (unsigned long be = 0; be < n_ent_elt; be++, ind_elt++) {

        if (version_fmt_gmsh < 41) {

          int  data[100];

          ecs_file_read(&data, sizeof(int),
                        1 + nbr_nod_elt_gmsh, fic_maillage);

          if (type_gmsh == GMSH_POINT1) {
            cpt_point += 1;
            continue;
          }
          else if (type_gmsh == GMSH_SEG2 || type_gmsh == GMSH_SEG3) {
            cpt_are += 1;
            continue;
          }

          label = data[0];

          /* Lecture des numéros des sommets de l'élément courant */

          for (ind_nod_elt = 0; ind_nod_elt < nbr_nod_elt_gmsh; ind_nod_elt++)
            num_nod_elt_gmsh[ind_nod_elt] = data[1 + ind_nod_elt];

        }
        else {

          size_t  data[100];

          ecs_file_read(&data, sizeof(size_t),
                        1 + nbr_nod_elt_gmsh, fic_maillage);

          if (type_gmsh == GMSH_POINT1) {
            cpt_point += 1;
            continue;
          }
          else if (type_gmsh == GMSH_SEG2 || type_gmsh == GMSH_SEG3) {
            cpt_are += 1;
            continue;
          }

          label = data[0];

          /* Lecture des numéros des sommets de l'élément courant */

          for (ind_nod_elt = 0; ind_nod_elt < nbr_nod_elt_gmsh; ind_nod_elt++)
            num_nod_elt_gmsh[ind_nod_elt] = data[1 + ind_nod_elt];

        }

        /* Stockage des valeurs avant transfert dans la structure `maillage' */
        /*===================================================================*/

        /* Identification de l'entité concernée */

        ent_num = ecs_maillage_pre__ret_typ_geo(type_ecs);

        /* Position des numéros de sommets du prochain élément */

        elt_pos_som_ent[ent_num][cpt_elt_ent[ent_num] + 1]
          = elt_pos_som_ent[ent_num][cpt_elt_ent[ent_num]] + nbr_som_elt;

        /* Connectivité de l'élément par ses numéros de sommets */

        for (ind_som_elt = 0; ind_som_elt < nbr_som_elt; ind_som_elt++) {

          elt_val_som_ent
            [ent_num]
            [elt_pos_som_ent[ent_num][cpt_elt_ent[ent_num]] - 1 + ind_som_elt]
            = num_nod_elt_gmsh
            [ecs_gmsh_elt_liste_c[type_gmsh - 1].num_som[ind_som_elt] - 1];

        }

        elt_val_fam_ent[ent_num][cpt_elt_ent[ent_num]]
          = tag_ent + fam_shift[ent_num];

        /* Incrémentation du nombre d'éléments lus */

        cpt_elt_ent[ent_num]++;

      }

      nbr_elt_lus += n_ent_elt;
    }

    ecs_file_set_type(fic_maillage, ECS_FILE_TYPE_TEXT);

    ecs_file_gets(chaine, 2, fic_maillage, num_ligne);

  }

  /* Indications */

  if (cpt_point > -1 || cpt_are > 0)
    printf("\n");
  if (cpt_point > -1)
    printf(_("  %ld elements of type \"point\" ignored\n"),
           (long)cpt_point);
  if (cpt_are > -1)
    printf(_("  %ld elements of type \"edge\" ignored\n"),
           (long)cpt_are);
  if (cpt_point > -1 || cpt_are > 0)
    printf("\n");

  /* Ligne de fin de rubrique */

  ecs_file_gets(chaine, ECS_LOC_LNG_MAX_CHAINE_GMSH, fic_maillage, num_ligne);

  if (   (strncmp(chaine, "$EndElements", strlen("$EndElements")) != 0)
      && (strncmp(chaine, "$ENDELM", strlen("$ENDELM")) != 0))
    ecs_error(__FILE__, __LINE__, 0,
              _("The end of the elements list specification of file\n"
                "\"%s\" is not present at the expected place.)"),
              ecs_file_get_name(fic_maillage));

  /* Réallocations des tableaux locaux */
  /*===================================*/

  /* Reallocations des tableaux locaux et mise à jour des références */
  /*=================================================================*/

  for (ient = ECS_ENTMAIL_FAC; ient < ECS_N_ENTMAIL; ient++) {

    if (cpt_elt_ent[ient] != 0) {
      ECS_REALLOC(elt_pos_som_ent[ient], cpt_elt_ent[ient] + 1, ecs_size_t);
      ECS_REALLOC(elt_val_som_ent[ient],
                  elt_pos_som_ent[ient][cpt_elt_ent[ient]] - 1, ecs_int_t);
      ECS_REALLOC(elt_val_fam_ent[ient], cpt_elt_ent[ient], int);
    }

    ecs_maillage_pre__label_en_indice
      (nbr_som,
       elt_pos_som_ent[ient][cpt_elt_ent[ient]] - 1,
       *som_val_label,
       elt_val_som_ent[ient]);
  }

  ECS_FREE(*som_val_label);

  /* Transfert des valeurs lues dans les structures d'entités de maillage */
  /*======================================================================*/

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

/*============================================================================
 *  Fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Lecture d'un fichier de maillage au format Gmsh
 *   et affectation des donnees dans la structure de maillage
 *----------------------------------------------------------------------------*/

ecs_maillage_t  *
ecs_pre_gmsh__lit_maillage(const char  *nom_fic_maillage)
{
  ecs_file_t   *fic_maillage;          /* Descripteur du fichier */
  char          chaine[ECS_LOC_LNG_MAX_CHAINE_GMSH]; /* Ligne lue  */
  int           num_ligne;             /* Compteur des lignes lues */
  int           version_fmt_gmsh;
  int           type_fmt_gmsh;
  ecs_int_t    *som_val_label;
  bool          bool_elements = false;
  bool          bool_noeuds = false;

  int           n_phys_names = 0;
  int           n_phys_names_max = 0;
  int          *phys_name_id = NULL;
  char        **phys_name = NULL;

  /* Création d'un maillage initialement vide (valeur de retour) */

  ecs_maillage_t  *maillage = ecs_maillage__cree_nodal();

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Affichage du titre */
  /*====================*/

  printf(_("\n\n"
           "Reading mesh from file in Gmsh format\n"
           "----------------------\n"));

  printf(_("  Mesh file: %s\n\n\n"),
         nom_fic_maillage);

  /* Initialisations */
  /*=================*/

  num_ligne = 0;

  som_val_label = NULL;

  int fam_shift[ECS_N_ENTMAIL] = {0, 0};

  /* Par défaut, la version est 1.0 en ASCII */

  version_fmt_gmsh = 10;
  type_fmt_gmsh = 0;

  /* Ouverture du fichier Gmsh en lecture */
  /*---------------------------------------*/

  fic_maillage = ecs_file_open(nom_fic_maillage,
                               ECS_FILE_MODE_READ,
                               ECS_FILE_TYPE_TEXT);

  /*================================================*/
  /* Boucle sur les lignes du fichier de maillage : */
  /* tant qu'on n'a pa atteint la fin du fichier    */
  /*================================================*/


  while (ecs_file_gets_try(chaine, ECS_LOC_LNG_MAX_CHAINE_GMSH,
                           fic_maillage, &num_ligne) != NULL) {

    /* Si la chaine lue est un descripteur de format */

    if (strncmp(chaine, "$MeshFormat", strlen("$MeshFormat")) == 0)

      ecs_loc_pre_gmsh__lit_vers_format(fic_maillage,
                                        &version_fmt_gmsh,
                                        &type_fmt_gmsh,
                                        &num_ligne);

    else if (strncmp(chaine,
                     "$PhysicalNames",
                     strlen("$PhysicalNames")) == 0)
      _read_physical_names(fic_maillage,
                           &num_ligne,
                           &n_phys_names,
                           &n_phys_names_max,
                           &phys_name_id,
                           &phys_name);

    else if (strncmp(chaine,
                     "$Entities",
                     strlen("$Entities")) == 0) {

      ecs_famille_t  *vect_famille[ECS_N_ENTMAIL];
      _read_entities(n_phys_names_max,
                     phys_name_id,
                     (const char**)phys_name,
                     fic_maillage,
                     &num_ligne,
                     version_fmt_gmsh,
                     type_fmt_gmsh,
                     fam_shift,
                     vect_famille);

      ecs_maillage__definit_famille(maillage, vect_famille);
    }

    /* Si la chaine lue marque le début de la section des noeuds */

    else if (strncmp(chaine,
                     "$ParametricNodes",
                     strlen("$ParametricNodes")) == 0)
      ecs_error
        (__FILE__, __LINE__, 0,
         _("Error reading file\n\"%s\" :\n"
           "Nodes are defined in the \"ParametricNodes\" variant, which\n"
           "is not handled (and not documented as of Gmsh 2.4.0)."),
         ecs_file_get_name(fic_maillage));

    /* Si la chaine lue marque le début de la section des noeuds */

    else if (   (strncmp(chaine, "$Nodes", strlen("$Nodes")) == 0)
             || (strncmp(chaine, "$NOD", strlen("$NOD")) == 0)) {

      if (version_fmt_gmsh/10 == 4)
        ecs_loc_pre_gmsh__lit_nodes_v4(maillage,
                                       version_fmt_gmsh,
                                       fic_maillage,
                                       type_fmt_gmsh,
                                       &num_ligne,
                                       &som_val_label);
      else
        ecs_loc_pre_gmsh__lit_nodes(maillage,
                                    fic_maillage,
                                    type_fmt_gmsh,
                                    &num_ligne,
                                    &som_val_label);

      bool_noeuds = true;

    }

    /* Si la chaine lue marque le début de la section des éléments */

    else if (   (strncmp(chaine, "$Elements", strlen("$Elements")) == 0)
             || (strncmp(chaine, "$ELM", strlen("$ELM")) == 0)) {

      if (version_fmt_gmsh/10 == 4)
        ecs_loc_pre_gmsh__lit_elements_v4(maillage,
                                          fic_maillage,
                                          version_fmt_gmsh,
                                          type_fmt_gmsh,
                                          &num_ligne,
                                          maillage->n_vertices,
                                          fam_shift,
                                          &som_val_label);
      else
        ecs_loc_pre_gmsh__lit_elements(maillage,
                                       fic_maillage,
                                       version_fmt_gmsh,
                                       type_fmt_gmsh,
                                       &num_ligne,
                                       maillage->n_vertices,
                                       &som_val_label);

      bool_elements = true;

    }
  }

  if (ecs_file_eof(fic_maillage) == 0)
    ecs_error(__FILE__, __LINE__, errno,
              _("Error reading line %ld of file \"%s\"."),
              (long)num_ligne, ecs_file_get_name(fic_maillage));

  /* else : la fin de fichier a bien ete atteinte */


  /* On verifie qu'on a bien lu des noeuds et des elements */

  if (bool_noeuds == false)
    ecs_error(__FILE__, __LINE__, 0,
              _("Error reading a Gmsh mesh file:\n"
                "at line %ld of file \"%s\".\n"
                "The points definition was not found."),
              (long)num_ligne, ecs_file_get_name(fic_maillage));

  if (bool_elements == false)
    ecs_error(__FILE__, __LINE__, 0,
              _("Error reading a Gmsh mesh file:\n"
                "at line %ld of file \"%s\".\n"
                "The elements definition was not found."),
              (long)num_ligne, ecs_file_get_name(fic_maillage));


  /* Fermeture du fichier de maillage */
  /*----------------------------------*/

  ecs_file_free(fic_maillage);

  /* Transformation des familles en attributs "groupe" */
  /*---------------------------------------------------*/

  if (version_fmt_gmsh/10 == 4)
    ecs_maillage__cree_attributs(maillage);

  for (int i = 0; i < n_phys_names; i++)
    ECS_FREE(phys_name[i]);

  ECS_FREE(phys_name);
  ECS_FREE(phys_name_id);

  /* Retour */
  /*========*/

  return maillage;
}

/*----------------------------------------------------------------------------*/

