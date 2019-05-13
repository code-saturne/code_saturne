/*============================================================================
 * Read mesh from STAR-CCM+ format file
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

#include "cs_config.h"

#if defined(HAVE_CCM)

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "ecs_def.h"
#include "ecs_elt_typ_liste.h"
#include "ecs_mem.h"
#include "ecs_tab.h"

#include "ecs_descr.h"
#include "ecs_descr_chaine.h"
#include "ecs_maillage.h"
#include "ecs_maillage_priv.h"

/*----------------------------------------------------------------------------
 *  Fichiers `include' visibles du  paquetage courant
 *----------------------------------------------------------------------------*/

#include "ecs_maillage_pre.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "ecs_pre_ccm.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif

#include <libccmio/ccmio.h>
#include <libccmio/ccmioversion.h>

#ifdef __cplusplus
}
#endif

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

#if (kCCMIOVersion == 20601)
typedef CCMIOSize CCMIOSize_t;
#endif

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/* Vertices definition */

typedef struct {
  size_t           n_nodes; /* Number of vertices */
  ecs_int_t       *id;      /* Vertex labels */
  ecs_coord_t     *coord;   /* Coordinates (interlaced) */
} _nodes_t;


/* Face definitions (boundary and internal) */

typedef struct {
  size_t        n_faces;        /* Number of */
  ecs_int_t    *nbr_n;          /* Number of vertices/face */
  size_t        taille_connect; /* Connectivity size */
  ecs_int_t    *connect;        /* Vertices connectivity */
  ecs_int_t    *icel1;          /* Cell on positive face side */
  ecs_int_t    *icel2;          /* Cell on negative face side */
  ecs_int_t    *icoul;          /* Face colors */
} _faces_t;

/* Cell definitions */

typedef struct {
  size_t           n_cells;    /* Number of cells */
  ecs_int_t       *id;         /* Cell ids */
  ecs_int_t       *icoul;      /* Cell colors (types) */
} _cells_t;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Read vertex data.
 *----------------------------------------------------------------------------*/

static void
_read_vertices(CCMIOID   vertices,
               _nodes_t *nodes)
{
  unsigned int i;
  float scale;
  double *verts_coo;
  CCMIOID map_id;
#if (kCCMIOVersion == 20601)
  int dims = 1;
  CCMIOSize n_vertices = 0, size;
#else
  CCMIOSize_t dims = 1, n_vertices = 0, size;
#endif
  int *map_data;
  CCMIOError err = kCCMIONoErr;

  /* Get description */

  CCMIOEntityDescription(&err, vertices, &size, NULL);

  if (err == kCCMIONoErr && size > 0) {
    char *s_descr;
    ECS_MALLOC(s_descr, size + 1, char);
    CCMIOEntityDescription(&err, vertices, &size, s_descr);
    if (err == kCCMIONoErr)
      printf(_("\n  Vertex coordinates set: %s\n"), s_descr);
    ECS_FREE(s_descr);
  }

  /* Read data */

  CCMIOEntitySize(&err, vertices, &n_vertices, NULL);
  if (err != kCCMIONoErr)
    ecs_error(__FILE__, __LINE__, 0,
              _("Error %d reading vertex information in STAR-CCM+ file."),
              (int)err);

  nodes->n_nodes = (ecs_int_t)n_vertices;

  ECS_MALLOC(verts_coo, 3*nodes->n_nodes, double);
  ECS_MALLOC(map_data, nodes->n_nodes, int);

  CCMIOReadVerticesd(&err, vertices, &dims, &scale, &map_id, verts_coo,
                     0, n_vertices);
  CCMIOReadMap(&err, map_id, map_data, 0, n_vertices);

  /* Transfer to local node structure */

  ECS_MALLOC(nodes->id, nodes->n_nodes, ecs_int_t);
  ECS_MALLOC(nodes->coord, 3*nodes->n_nodes, ecs_coord_t);

  for (i = 0; i < n_vertices; i++) {

    nodes->id[i] = map_data[i];

    nodes->coord[dims*i  ] = scale*verts_coo[dims*i];
    nodes->coord[dims*i+1] = scale*verts_coo[dims*i+1];
    nodes->coord[dims*i+2] = scale*verts_coo[dims*i+2];

  }

  ECS_FREE(map_data);
  ECS_FREE(verts_coo);
}

/*----------------------------------------------------------------------------
 * Read cell data.
 *----------------------------------------------------------------------------*/

static void
_read_cells(CCMIOID   topology,
            _cells_t *cells)
{
  unsigned int i;

  CCMIOID map_id, id;
  CCMIOSize_t n_cells;
  CCMIOError err = kCCMIONoErr;

  int *map_data;
  int *cell_type;

  static const int k_cell_inc = 4;

  CCMIOGetEntity(&err, topology, kCCMIOCells, 0, &id);
  CCMIOEntitySize(&err, id, &n_cells, NULL);

  ECS_MALLOC(cell_type, n_cells, int);
  ECS_MALLOC(map_data, n_cells, int);

  for (i = 0;  i < n_cells;  i += k_cell_inc) {
    CCMIOReadCells(&err, id, &map_id, &cell_type[i], i, i + k_cell_inc);
    CCMIOReadMap(&err, map_id, &map_data[i], i, i + k_cell_inc);
  }

  /* Transfer to local cell structure */

  cells->n_cells = (ecs_int_t) n_cells;

  ECS_MALLOC(cells->id, n_cells, ecs_int_t);
  ECS_MALLOC(cells->icoul, n_cells, ecs_int_t);

  for (i = 0; i < n_cells; i++) {
    cells->id[i] = map_data[i];
    cells->icoul[i] = cell_type[i];
  }

  ECS_FREE(map_data);
  ECS_FREE(cell_type);
}

/*----------------------------------------------------------------------------
 * Read mesh.
 *----------------------------------------------------------------------------*/

static void
_read_mesh(CCMIOID   vertices,
           CCMIOID   topology,
           _nodes_t *nodes,
           _faces_t *i_faces,
           _faces_t *b_faces,
           _cells_t *cells)
{
  size_t i;
  int j;
  CCMIOID map_id, id;
  CCMIOError err = kCCMIONoErr;
  CCMIOSize_t n_faces, size;
  int *face_nverts; int *face_cells;
  size_t pos, cpt;

  _read_vertices(vertices,
                 nodes);

  _read_cells(topology,
              cells);

  /* Read the topology description. */

  CCMIOEntityDescription(&err, topology, &size, NULL);

  if (err == kCCMIONoErr && size > 0) {
    char *s_descr;
    ECS_MALLOC(s_descr, size + 1, char);
    CCMIOEntityDescription(&err, topology, &size, s_descr);
    if (err == kCCMIONoErr)
      printf(_("  Face-based topology: %s\n"), s_descr);
    ECS_FREE(s_descr);
  }

  /* Read the internal faces. */

  CCMIOGetEntity(&err, topology, kCCMIOInternalFaces, 0, &id);

  CCMIOEntityDescription(&err, id, &size, NULL);

  if (err == kCCMIONoErr && size > 0) {
    char *s_descr;
    ECS_MALLOC(s_descr, size + 1, char);
    CCMIOEntityDescription(&err, id, &size, s_descr);
    if (err == kCCMIONoErr)
      printf(_("    Faces: %s\n"), s_descr);
    ECS_FREE(s_descr);
  }

  /* Read internal faces data */

  CCMIOEntitySize(&err, id, &n_faces, NULL);

  ECS_MALLOC(face_cells, 2*n_faces, int);

  CCMIOReadFaces(&err, id, kCCMIOInternalFaces, NULL, &size, NULL,
                 kCCMIOStart, kCCMIOEnd);

  if (err != kCCMIONoErr)
    ecs_error(__FILE__, __LINE__, 0,
              _("Error %d reading faces information in STAR-CCM+ file."),
              (int)err);

  ECS_MALLOC(face_nverts, size, int);

  CCMIOReadFaces(&err, id, kCCMIOInternalFaces, &map_id, NULL, face_nverts,
                 kCCMIOStart, kCCMIOEnd);

  if (err != kCCMIONoErr)
    ecs_error(__FILE__, __LINE__, 0,
              _("Error %d reading internal faces in STAR-CCM+ file."),
              (int)err);

  CCMIOReadFaceCells(&err, id, kCCMIOInternalFaces, face_cells,
                     kCCMIOStart, kCCMIOEnd);

  if (err != kCCMIONoErr)
    ecs_error(__FILE__, __LINE__, 0,
              _("Error %d reading face->cell connectivity in STAR-CCM+ file."),
              (int)err);

  /* Map data is not used here; if it were necessary, we would have:
   *
   * int *map_data;
   *
   * ECS_MALLOC(map_data, n_faces, int);
   * CCMIOReadMap(&err, map_id, map_data, kCCMIOStart, kCCMIOEnd);
   * ECS_FREE(map_data);
   */

  if (err != kCCMIONoErr)
    ecs_error(__FILE__, __LINE__, 0,
              _("Error %d reading map in STAR-CCM+ file."),
              (int)err);

  /* Transfer to local face structure */

  i_faces->n_faces = n_faces;
  i_faces->taille_connect = size;

  ECS_MALLOC(i_faces->nbr_n, n_faces, ecs_int_t);

  ECS_MALLOC(i_faces->icel1, n_faces, ecs_int_t);
  ECS_MALLOC(i_faces->icel2, n_faces, ecs_int_t);
  ECS_MALLOC(i_faces->icoul, n_faces, ecs_int_t);

  i_faces->taille_connect = size - n_faces;

  ECS_MALLOC(i_faces->connect, i_faces->taille_connect, ecs_int_t);

  for (i = 0, pos = 0, cpt = 0; i < i_faces->n_faces; i++) {

    i_faces->nbr_n[i] = face_nverts[pos];
    i_faces->icel1[i] = face_cells[2*i];
    i_faces->icel2[i] = face_cells[2*i+1];
    assert(i_faces->icel1[i] > 0);
    assert(i_faces->icel2[i] > 0);
    i_faces->icoul[i] = 0;

    for (j = 0; j < i_faces->nbr_n[i]; j++)
      i_faces->connect[cpt + j] = (ecs_int_t) face_nverts[pos + j + 1];

    cpt += i_faces->nbr_n[i];
    pos += face_nverts[pos] + 1;
  }

  ECS_FREE(face_nverts);
  ECS_FREE(face_cells);

  /* Read the boundary faces */

  {
#if (kCCMIOVersion == 20601)
    int _index = 0;
#else
    CCMIOSize_t _index = 0;
#endif
    int i_beg = 0;
    int i_beg_c = 0;

    b_faces->n_faces = 0;
    b_faces->taille_connect = 0;

    while (   CCMIONextEntity(NULL, topology, kCCMIOBoundaryFaces, &_index, &id)
           == kCCMIONoErr) {

      int boundary_id;

      CCMIOEntityDescription(&err, id, &size, NULL);

      if (err == kCCMIONoErr && size > 0) {
        char *s_descr;
        ECS_MALLOC(s_descr, size + 1, char);
        CCMIOEntityDescription(&err, id, &size, s_descr);
        if (err == kCCMIONoErr)
          printf(_("    Faces: %s\n"), s_descr);
        ECS_FREE(s_descr);
      }

      CCMIOEntitySize(&err, id, &n_faces, NULL);

      ECS_REALLOC(face_cells, 2*n_faces, int);

      CCMIOReadFaces(&err, id, kCCMIOBoundaryFaces, NULL, &size, NULL,
                     kCCMIOStart, kCCMIOEnd);

      ECS_MALLOC(face_nverts, size, int);

      CCMIOReadFaces(&err, id, kCCMIOBoundaryFaces, &map_id, NULL,
                     face_nverts, kCCMIOStart, kCCMIOEnd);
      CCMIOReadFaceCells(&err, id, kCCMIOBoundaryFaces, face_cells,
                         kCCMIOStart, kCCMIOEnd);

      /* Map data is not used here; if it were necessary, we would have:
       *
       * int *map_data;
       *
       * ECS_MALLOC(map_data, n_faces, int);
       * CCMIOReadMap(&err, map_id, map_data, kCCMIOStart, kCCMIOEnd);
       * ECS_FREE(map_data);
      */

      CCMIOGetEntityIndex(&err, id, &boundary_id);

      /* Transfer to local face structure */

      b_faces->n_faces += n_faces;
      b_faces->taille_connect += (size - n_faces);

      if (i_beg == 0) {
        ECS_MALLOC(b_faces->nbr_n, b_faces->n_faces, ecs_int_t);
        ECS_MALLOC(b_faces->icel1, b_faces->n_faces, ecs_int_t);
        ECS_MALLOC(b_faces->icel2, b_faces->n_faces, ecs_int_t);
        ECS_MALLOC(b_faces->connect, b_faces->taille_connect, ecs_int_t);
        ECS_MALLOC(b_faces->icoul, b_faces->n_faces, ecs_int_t);
      }
      else {
        ECS_REALLOC(b_faces->nbr_n, b_faces->n_faces, ecs_int_t);
        ECS_REALLOC(b_faces->icel1, b_faces->n_faces, ecs_int_t);
        ECS_REALLOC(b_faces->icel2, b_faces->n_faces, ecs_int_t);
        ECS_REALLOC(b_faces->connect, b_faces->taille_connect, ecs_int_t);
        ECS_REALLOC(b_faces->icoul, b_faces->n_faces, ecs_int_t);
      }

      for (i = 0, pos = 0, cpt = 0; i < (size_t)n_faces; i++) {

        b_faces->nbr_n[i_beg + i] = face_nverts[pos];
        b_faces->icel1[i_beg + i] = face_cells[i];
        b_faces->icel2[i_beg + i] = -1;
        assert(b_faces->icel1[i_beg + i] > 0);
        b_faces->icoul[i_beg + i] = boundary_id;

        for (j = 0; j < b_faces->nbr_n[i_beg + i]; j++)
          b_faces->connect[i_beg_c + cpt + j] = face_nverts[pos + j + 1];

        cpt += b_faces->nbr_n[i_beg + i];
        pos += face_nverts[pos] + 1;
      }

      i_beg += n_faces;
      i_beg_c += (size - n_faces);

      ECS_FREE(face_nverts);
      ECS_FREE(face_cells);
    }
  }

  printf("\n");
}

/*----------------------------------------------------------------------------
 * Extract a mesh's cell -> faces connectivity.
 *
 * We consider a common numbering for interior and boundary faces, in which
 * boundary faces are defined first. The common id of the i-th boundary
 * face is thus i, while that of the j-th interior face is nbr_fbr + j.
 *----------------------------------------------------------------------------*/

static void
ecs_pre_ccm__cel_fac(const _faces_t         *faces,
                     const _cells_t         *cels,
                     ecs_int_t       **const p_pos_cel_fac,
                     ecs_int_t       **const p_val_cel_fac)
{
  ecs_int_t    icel, icel1, icel2, n_cells_loc;
  size_t       ifac;

  ecs_int_t  * cpt_cel_fac = NULL;
  ecs_int_t  * pos_cel_fac = NULL;
  ecs_int_t  * val_cel_fac = NULL;

  /* Allocate and initialize position index */

  n_cells_loc = cels->n_cells;

  ECS_MALLOC(pos_cel_fac, n_cells_loc + 1, ecs_int_t);

  for (icel = 0; icel < n_cells_loc + 1; icel++)
    pos_cel_fac[icel] = 0;

  /* Count number of faces per cell
   * (we assign the temporary counter for icel to pos_cel_fac[icel + 1]
   * instead of pos_cel_fac[icel] so as to simplify the next step) */

  for (ifac = 0; ifac < faces->n_faces; ifac++) {
    icel1 = faces->icel1[ifac] - 1;
    icel2 = faces->icel2[ifac] - 1;
    pos_cel_fac[icel1 + 1] += 1;
    if (icel2 >= 0)
      pos_cel_fac[icel2 + 1] += 1;
  }

  /* Build position index */

  pos_cel_fac[0] = 1;
  for (icel = 0; icel < n_cells_loc; icel++)
    pos_cel_fac[icel + 1] = pos_cel_fac[icel] + pos_cel_fac[icel + 1];

  /* Build array of values */

  ECS_MALLOC(val_cel_fac, pos_cel_fac[n_cells_loc] - 1, ecs_int_t);
  ECS_MALLOC(cpt_cel_fac, n_cells_loc, ecs_int_t);

  for (icel = 0; icel < n_cells_loc; icel++)
    cpt_cel_fac[icel] = 0;

  for (ifac = 0; ifac < faces->n_faces; ifac++) {
    icel1 = faces->icel1[ifac] - 1;
    icel2 = faces->icel2[ifac] - 1;
    assert(icel1 >= 0);
    val_cel_fac[pos_cel_fac[icel1] + cpt_cel_fac[icel1] - 1] = (ifac + 1);
    cpt_cel_fac[icel1] += 1;
    if (icel2 >= 0) {
      val_cel_fac[pos_cel_fac[icel2] + cpt_cel_fac[icel2] - 1] = -(ifac + 1);
      cpt_cel_fac[icel2] += 1;
    }
  }

  ECS_FREE(cpt_cel_fac);

  /* Return values */

  *p_pos_cel_fac = pos_cel_fac;
  *p_val_cel_fac = val_cel_fac;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  {
    ecs_int_t ipos, ival;

    printf("dbg : cs_loc_maillage_ret_cel_fac\n"
           "nombre de cellules extraites = %d\n", n_cells_extr);
    for (ipos = 0; ipos < n_cells_extr; ipos++) {
      printf("  cellule %d\n", ipos);
      printf("    pos_cel_fac[%d] = %d\n", ipos, pos_cel_fac[ipos]);
      for (ival = pos_cel_fac[ipos]     - 1;
           ival < pos_cel_fac[ipos + 1] - 1;
           ival++)
        printf("      val_cel_fac[%d] = %d\n", ival, val_cel_fac[ival]);
    }
    printf("  pos_cel_fac[%d] = %d\n", ipos, pos_cel_fac[ipos]);
  }
#endif

}

/*----------------------------------------------------------------------------
 * Build mesh structure
 *----------------------------------------------------------------------------*/

static ecs_maillage_t *
ecs_pre_ccm__prepa_mail(_nodes_t    *noeuds,
                        _faces_t    *faces,
                        _cells_t    *cels,
                        ecs_int_t   *pos_cel_fac,
                        ecs_int_t   *val_cel_fac)
{
  /* Temporary arrays used before transfer to mesh structure */

  size_t       cpt_elt_ent        [ECS_N_ENTMAIL]; /* N. elts per entity */
  ecs_int_t    cpt_coul_ent       [ECS_N_ENTMAIL]; /* Color counter */
  size_t       cpt_val_som_ent    [ECS_N_ENTMAIL]; /* Connectivity size */
  ecs_int_t  * val_coul_ent       [ECS_N_ENTMAIL]; /* Array of colors */
  ecs_size_t * cpt_elt_coul_ent   [ECS_N_ENTMAIL];
  ecs_size_t * elt_pos_som_ent    [ECS_N_ENTMAIL]; /* Vertex positions */
  ecs_int_t  * elt_val_som_ent    [ECS_N_ENTMAIL]; /* Vertex numbers */
  ecs_int_t  * elt_val_color_ent  [ECS_N_ENTMAIL]; /* Element colors */

  size_t ifac;
  ecs_int_t nbr_som_elt;
  ecs_int_t icoul;

  ecs_entmail_t    entmail_e;

  ecs_int_t      ient;
  size_t         pos_elt;
  ecs_int_t      isom;
  size_t         n_faces;
  size_t         icel;
  size_t         pos_fac;
  ecs_int_t      nbr_som_fac;
  ecs_int_t      num_som_deb_fac;

  /* Create intially empty mesh (return value) */

  ecs_maillage_t  *maillage = ecs_maillage__cree_nodal();

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Vertices */

  ecs_maillage_pre__cree_som(maillage,
                             noeuds->n_nodes,
                             noeuds->coord);

  /* Cells */

  for (ient = 0; ient < ECS_N_ENTMAIL; ient++) {

    cpt_elt_ent        [ient] = 0;
    cpt_val_som_ent    [ient] = 0;

    elt_pos_som_ent    [ient] = NULL;
    elt_val_som_ent    [ient] = NULL;
    elt_val_color_ent  [ient] = NULL;

    cpt_coul_ent       [ient] = 0;
    val_coul_ent       [ient] = NULL;
    cpt_elt_coul_ent   [ient] = NULL;

  }

  /* Faces (polygons) */

  entmail_e = ECS_ENTMAIL_FAC;

  cpt_elt_ent[entmail_e] = faces->n_faces;
  cpt_val_som_ent[entmail_e] = faces->taille_connect;
  elt_val_som_ent[entmail_e] = faces->connect;

  if (cpt_elt_ent[entmail_e] > 0) {

    ECS_MALLOC(elt_pos_som_ent[entmail_e],
               cpt_elt_ent[entmail_e] + 1,
               ecs_size_t);

    elt_pos_som_ent[entmail_e][0] = 1;

    ECS_MALLOC(elt_val_color_ent[entmail_e],
               cpt_elt_ent[entmail_e],
               ecs_int_t);
  }

  for (ifac = 0; ifac < faces->n_faces; ifac++) {

    /* Test if color is assigned */

    for (icoul = 0;
            icoul < cpt_coul_ent[entmail_e]
         && val_coul_ent[entmail_e][icoul] != faces->icoul[ifac];
         icoul++);

    if (icoul == cpt_coul_ent[entmail_e]) {

      /* This color value has not been saved yet */

      ECS_REALLOC(val_coul_ent[entmail_e],
                  cpt_coul_ent[entmail_e] + 1,
                  ecs_int_t);
      ECS_REALLOC(cpt_elt_coul_ent[entmail_e],
                  cpt_coul_ent[entmail_e] + 1,
                  ecs_size_t);

      cpt_elt_coul_ent[entmail_e][icoul] = 0;
      val_coul_ent[entmail_e][icoul] = faces->icoul[ifac];
      cpt_coul_ent[entmail_e]++;

    }

    /* Assign color */

    cpt_elt_coul_ent[entmail_e][icoul]++;
    elt_val_color_ent[entmail_e][ifac] = icoul + 1;

    nbr_som_elt = faces->nbr_n[ifac];

    /* Build connectivity */

    pos_elt = elt_pos_som_ent[entmail_e][ifac];

    elt_pos_som_ent[entmail_e][ifac + 1] = pos_elt + nbr_som_elt;

  }

  /* Cells (polyhedra) */

  entmail_e = ECS_ENTMAIL_CEL;

  cpt_elt_ent[entmail_e] = cels->n_cells;

  n_faces = pos_cel_fac[cels->n_cells] - 1;

  cpt_val_som_ent[entmail_e] = 0;

  for (ifac = 0; ifac < n_faces; ifac++)
    cpt_val_som_ent[entmail_e] += faces->nbr_n[ECS_ABS(val_cel_fac[ifac])-1] + 1;

  if (cpt_elt_ent[entmail_e] > 0) {

    ECS_MALLOC(elt_pos_som_ent[entmail_e],
               cpt_elt_ent[entmail_e] + 1,
               ecs_size_t);

    elt_pos_som_ent[entmail_e][0] = 1;

    ECS_MALLOC(elt_val_som_ent[entmail_e],
               cpt_val_som_ent[entmail_e],
               ecs_int_t);

    ECS_MALLOC(elt_val_color_ent[entmail_e],
               cpt_elt_ent[entmail_e],
               ecs_int_t);

  }

  for (icel = 0; icel < cels->n_cells; icel++){

    /* Position in connectivity */

    cpt_val_som_ent[entmail_e] = 0;

    for (ifac = pos_cel_fac[icel] - 1;
         ifac < (size_t)(pos_cel_fac[icel + 1] - 1);
         ifac++)
      cpt_val_som_ent[entmail_e]
        += faces->nbr_n[ECS_ABS(val_cel_fac[ifac]) - 1] + 1;

    elt_pos_som_ent[entmail_e][icel + 1]
      = elt_pos_som_ent[entmail_e][icel] + cpt_val_som_ent[entmail_e];

    /* Test if color is assigned */

    for (icoul = 0;
            icoul < cpt_coul_ent[entmail_e]
         && val_coul_ent[entmail_e][icoul] != cels->icoul[icel];
         icoul++);

    if (icoul == cpt_coul_ent[entmail_e]) {

      /* This color value has not been saved yet */

      ECS_REALLOC(val_coul_ent[entmail_e],
                  cpt_coul_ent[entmail_e] + 1,
                  ecs_int_t);
      ECS_REALLOC(cpt_elt_coul_ent[entmail_e],
                  cpt_coul_ent[entmail_e] + 1,
                  ecs_size_t);

      cpt_elt_coul_ent[entmail_e][icoul] = 0;
      val_coul_ent[entmail_e][icoul] = cels->icoul[icel];
      cpt_coul_ent[entmail_e]++;

    }

    /* Assign color */

    cpt_elt_coul_ent[entmail_e][icoul]++;
    elt_val_color_ent[entmail_e][icel] = icoul + 1;

  }

  /* Connectivity */

  /* Loop on all faces of val_cel_fac */

  cpt_val_som_ent[entmail_e] = 0;

  for (ifac = 0; ifac < n_faces; ifac++) {

    pos_fac = elt_pos_som_ent[ECS_ENTMAIL_FAC][ECS_ABS(val_cel_fac[ifac])-1];

    nbr_som_fac = faces->nbr_n[ECS_ABS(val_cel_fac[ifac])-1];

    /* Orientation */

    if (val_cel_fac[ifac] < 0) {

      num_som_deb_fac = faces->connect[pos_fac + (nbr_som_fac - 1) -1];

      for (isom = 0; isom < nbr_som_fac; isom++)
        elt_val_som_ent[entmail_e][cpt_val_som_ent[entmail_e] + isom]
          = faces->connect[pos_fac -1 + ((nbr_som_fac - 1) - isom)];

    }

    else {

      num_som_deb_fac = faces->connect[pos_fac - 1];

      for (isom = 0; isom < nbr_som_fac; isom++)
        elt_val_som_ent[entmail_e][cpt_val_som_ent[entmail_e] + isom]
          = faces->connect[pos_fac -1 + isom];

    }

    elt_val_som_ent[entmail_e][cpt_val_som_ent[entmail_e] + isom]
      = num_som_deb_fac;

    cpt_val_som_ent[entmail_e] += nbr_som_fac + 1;

  }

  /* Transfer values read to mesh entity structures */

  ecs_maillage_pre__cree_elt(maillage,
                             cpt_elt_ent,
                             elt_pos_som_ent,
                             elt_val_som_ent,
                             NULL,
                             elt_val_color_ent,
                             cpt_coul_ent,
                             val_coul_ent,
                             cpt_elt_coul_ent);

  ECS_FREE(pos_cel_fac);
  ECS_FREE(val_cel_fac);

  return maillage;
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Read of mesh from a STAR-CCM+ format file.
 *
 * parameters:
 *   nom_fic_maillage <-- mesh file name
 *   num_maillage     <-- mesh number
 *----------------------------------------------------------------------------*/

ecs_maillage_t *
ecs_pre_ccm__lit_maillage(const char  *nom_fic_maillage,
                          int          num_maillage)
{
  ecs_maillage_t * maillage = NULL;

  CCMIOID root, state, processor, vertices, topology;
  CCMIOError err;

#if (kCCMIOVersion == 20601)
  int i = 0;
  unsigned int size;
#else
  CCMIOSize_t i = 0;
  CCMIOSize_t size = 0;
#endif

  _nodes_t *nodes;
  _faces_t *liste_faces_bord;
  _faces_t *liste_faces_internes;
  _faces_t *liste_faces;
  _cells_t *liste_cels;

  ecs_tab_int_t tab_connect;
  ecs_tab_int_t tab_label;

  size_t    nbr_fabord;
  size_t    nbr_faint;
  size_t    n_faces;
  size_t    isom, ifac;

  ecs_int_t *pos_cel_fac;
  ecs_int_t *val_cel_fac;

  int       n_states = 0;
  size_t    taille_connect_bord = 0;
  size_t    taille_connect_interne = 0;
  size_t    taille_connect = 0;

  int   version = 0, maj_ver = 0, min_ver = 0, rel_ver = 0;

  CCMIOID *states = NULL;
  char *title = NULL;

  /*Xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  err = CCMIOOpenFile(NULL, nom_fic_maillage, kCCMIORead, &root);

  if (err != kCCMIONoErr)
    ecs_error(__FILE__, __LINE__, 0,
              _("Error %d opening STAR-CCM+ file: %s"),
              (int)err, nom_fic_maillage);

  CCMIOGetVersion(&err, root.node, &version);

  maj_ver =  version / 10000;
  min_ver = (version % 10000)/ 100;
  rel_ver = (version % 100);

  CCMIOGetTitle(&err, root.node, &title);

  printf(_("\n\n"
           "Reading mesh from file in STAR-CCM+ format\n"
           "----------------------\n"));

  printf(_("  Mesh file: %s\n\n\n"
           "  STAR-CCM+ format version: %1d.%1d.%1d\n"
           "  Title: %s\n\n"),
             nom_fic_maillage, maj_ver, min_ver, rel_ver, title);

  /* The string 'title' has been allocated inside CCMIOGetTitle */

  free(title);

  /* Initialisations */

  nodes = NULL;
  liste_faces_bord  = NULL;
  liste_faces_internes  = NULL;
  liste_cels  = NULL;

  ECS_MALLOC(nodes, 1, _nodes_t);
  ECS_MALLOC(liste_faces_bord, 1, _faces_t);
  ECS_MALLOC(liste_faces_internes, 1, _faces_t);
  ECS_MALLOC(liste_cels, 1, _cells_t);

  /* Search for available states */

  err = kCCMIONoErr;

  printf(_("  Available meshes:\n"));

  while (err == kCCMIONoErr) {

    char name[kCCMIOMaxStringLength + 1];

    CCMIONextEntity(&err, root, kCCMIOState, &i, &state);

    if (err != kCCMIONoErr)
      break;

    CCMIOEntityName(&err, state, name);

    if (err != kCCMIONoErr)
      ecs_error(__FILE__, __LINE__, 0,
                _("Error %d reading state name in STAR-CCM+ file."),
                (int)err);

    printf("    %2d: \"%s\"\n", n_states + 1, name);

    ECS_REALLOC(states, n_states + 1, CCMIOID);
    states[n_states] = state;

    n_states += 1;
  }

  err = kCCMIONoErr;

  if (num_maillage == 0) {
    printf(_("\n  No mesh was specified; the first is read.\n\n"));
    state = states[0];
  }
  else {
    printf(_("\n  Mesh number %d was specified.\n\n"), num_maillage);
    if (num_maillage > 0 && num_maillage <= n_states)
      state = states[num_maillage - 1];
    else if (num_maillage > 0)
      ecs_error(__FILE__, __LINE__, 0,
              _("The specified mesh number (%d) is greater than\n"
                "the number of meshes defined (%d) in file\n%s.\n"),
              num_maillage, n_states, nom_fic_maillage);
    else
      ecs_error(__FILE__, __LINE__, 0,
              _("The specified mesh number (%d) is negative."),
              num_maillage);
  }

  ECS_FREE(states);

  CCMIOEntityDescription(&err, state, &size, NULL);

  if (err != kCCMIONoErr)
    ecs_error(__FILE__, __LINE__, 0,
              _("Error %d reading entity description in STAR-CCM+ file."),
              (int)err);

  if (size > 0) {
    char *desc;
    ECS_MALLOC(desc, size + 1, char);
    CCMIOEntityDescription(&err, state, NULL, desc);
    printf("  Mesh/state description: %s)\n", desc);
    ECS_FREE(desc);
  }

  i = 0;
  CCMIONextEntity(&err, state, kCCMIOProcessor, &i, &processor);

  if (err != kCCMIONoErr)
    ecs_error(__FILE__, __LINE__, 0,
              _("Error %d reading processor information in STAR-CCM+ file."),
              (int)err);

  CCMIOReadProcessor(&err, processor, &vertices, &topology, NULL, NULL);

  if (err != kCCMIONoErr)
    ecs_error(__FILE__, __LINE__, 0,
              _("Error %d reading vertices and topology information\n"
                "in STAR-CCM+ file."), (int)err);

  _read_mesh(vertices,
             topology,
             nodes,
             liste_faces_internes,
             liste_faces_bord,
             liste_cels);

  /* Fermetures */

  CCMIOCloseFile(&err, vertices);
  CCMIOCloseFile(&err, topology);

  CCMIOCloseFile(&err, root);

  /* Switch to preprocessor type structure */

  /* switching to indexes is done before appending faces so as to keep
     icel2 = -1 for boundary faces */

  /* vertex labels */

  /* boundary faces */

  tab_connect.nbr = liste_faces_bord->taille_connect;
  tab_connect.val = liste_faces_bord->connect;

  tab_label.nbr = nodes->n_nodes;
  tab_label.val = nodes->id;

  ecs_tab_int__ref_en_indice(tab_connect, tab_label, false);
  liste_faces_bord->connect = tab_connect.val;

  /* interior faces */

  tab_connect.nbr = liste_faces_internes->taille_connect;
  tab_connect.val = liste_faces_internes->connect;

  ecs_tab_int__ref_en_indice(tab_connect, tab_label, false);
  liste_faces_internes->connect = tab_connect.val;

  /* cell labels */

  tab_label.nbr = liste_cels->n_cells;
  tab_label.val = liste_cels->id;
  tab_connect.nbr = liste_faces_bord->n_faces;

  /* For boundary faces, only icel1 is handled (icel2 = -1) */

  tab_connect.val = liste_faces_bord->icel1;
  ecs_tab_int__ref_en_indice(tab_connect, tab_label, false);
  liste_faces_bord->icel1 = tab_connect.val;

  /* Interior faces */

  tab_connect.nbr = liste_faces_internes->n_faces;

  tab_connect.val = liste_faces_internes->icel1;
  ecs_tab_int__ref_en_indice(tab_connect, tab_label, false);
  liste_faces_internes->icel1 = tab_connect.val;

  tab_connect.val = liste_faces_internes->icel2;
  ecs_tab_int__ref_en_indice(tab_connect, tab_label, false);
  liste_faces_internes->icel2 = tab_connect.val;

  /* shift ... n-1 to 1 ... n */

  for (isom = 0; isom < liste_faces_internes->taille_connect; isom++)
    liste_faces_internes->connect[isom] += 1;

  for (isom = 0; isom < liste_faces_bord->taille_connect; isom++)
    liste_faces_bord->connect[isom] += 1;

  for (ifac = 0; ifac < liste_faces_internes->n_faces; ifac++){
    assert(liste_faces_internes->icel1[ifac] > -1);
    liste_faces_internes->icel1[ifac] += 1;
    if (liste_faces_internes->icel2[ifac] != -1)
      liste_faces_internes->icel2[ifac] += 1;
  }

  for (ifac = 0; ifac < liste_faces_bord->n_faces; ifac++)
    liste_faces_bord->icel1[ifac] += 1;

  /* Append boundary and interior faces, with interior faces first */

  taille_connect_bord = liste_faces_bord->taille_connect;
  taille_connect_interne = liste_faces_internes->taille_connect;

  taille_connect = taille_connect_bord + taille_connect_interne;

  nbr_fabord = liste_faces_bord->n_faces;
  nbr_faint  = liste_faces_internes->n_faces;

  n_faces = nbr_fabord + nbr_faint;

  ECS_REALLOC(liste_faces_bord->connect, taille_connect, ecs_int_t);
  ECS_REALLOC(liste_faces_bord->icoul, n_faces, ecs_int_t);
  ECS_REALLOC(liste_faces_bord->icel1, n_faces, ecs_int_t);
  ECS_REALLOC(liste_faces_bord->icel2, n_faces, ecs_int_t);
  ECS_REALLOC(liste_faces_bord->nbr_n, n_faces, ecs_int_t);

  liste_faces = liste_faces_bord;

  for (isom = 0; isom < taille_connect_interne; isom++)
    liste_faces->connect[taille_connect_bord + isom ]
      = liste_faces_internes->connect[isom];

  for (ifac = 0; ifac < nbr_faint; ifac ++){
    liste_faces->icoul[nbr_fabord + ifac] = liste_faces_internes->icoul[ifac];
    liste_faces->icel1[nbr_fabord + ifac] = liste_faces_internes->icel1[ifac];
    liste_faces->icel2[nbr_fabord + ifac] = liste_faces_internes->icel2[ifac];
    liste_faces->nbr_n[nbr_fabord + ifac] = liste_faces_internes->nbr_n[ifac];
  }

  liste_faces->taille_connect = taille_connect;
  liste_faces->n_faces = n_faces;

  /* Libération de la memoire */

  ECS_FREE(liste_faces_internes->icoul);
  ECS_FREE(liste_faces_internes->icel1);
  ECS_FREE(liste_faces_internes->icel2);
  ECS_FREE(liste_faces_internes->nbr_n);
  ECS_FREE(liste_faces_internes->connect);

  ECS_FREE(liste_faces_internes);

  /* Préparation de la connectivité cellule */

  ecs_pre_ccm__cel_fac(liste_faces,
                       liste_cels,
                       &pos_cel_fac,
                       &val_cel_fac);

  /* Remplissage des tableaux */

  maillage = ecs_pre_ccm__prepa_mail(nodes,
                                     liste_faces,
                                     liste_cels,
                                     pos_cel_fac,
                                     val_cel_fac);

  ECS_FREE(nodes->id);
  ECS_FREE(nodes);

  ECS_FREE(liste_faces->icel1);
  ECS_FREE(liste_faces->icel2);
  ECS_FREE(liste_faces->nbr_n);
  ECS_FREE(liste_faces->icoul);

  ECS_FREE(liste_faces);

  ECS_FREE(liste_cels->icoul);
  ECS_FREE(liste_cels->id);

  ECS_FREE(liste_cels);

  /* Renvoi de la structure de maillage */

  return maillage;
}

/*----------------------------------------------------------------------------*/

#endif /* HAVE_CCM */
