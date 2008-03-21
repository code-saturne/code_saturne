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
 * Manage the exchange of data between Code_Saturne and the pre-processor
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_error.h>
#include <bft_mem.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <fvm_periodicity.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_mesh.h"
#include "cs_perio.h"
#include "cs_pp_io.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_ecs_messages.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Add a periodicity to mesh->periodicities (fvm_periodicity_t *) structure.
 * Make the call fvm_periodicity_add_by_matrix()
 *
 * Returns:
 *----------------------------------------------------------------------------*/

static void _add_periodicity(void)
{
  cs_int_t  i, j, tr_id;
  double  matrix[3][4];

  cs_mesh_builder_t  *mesh_builder = cs_glob_mesh_builder;
  cs_mesh_t  *mesh = cs_glob_mesh;

  cs_bool_t  inv_pt = CS_FALSE;

  assert(sizeof(cs_real_t) == sizeof(double));

  for (i = 0; i < CS_DIM_3; i++)
    if (CS_ABS(mesh_builder->invariant_point[i]) >= 1.e-12)
      inv_pt = CS_TRUE;

  if (mesh_builder->perio_type == FVM_PERIODICITY_TRANSLATION) {

    bft_printf(_(" Ajout d'une périodicité en translation  ..."));
    bft_printf_flush();

    tr_id = fvm_periodicity_add_translation(mesh->periodicity,
                                            mesh_builder->perio_num,
                                            mesh_builder->translation);

    bft_printf("\t\t [ok]\n");

  }
  else if (mesh_builder->perio_type == FVM_PERIODICITY_ROTATION) {

    bft_printf(_(" Ajout d'une périodicité en rotation  ..."));
    bft_printf_flush();

    for (i = 0; i < 3; i++)
      for (j = 0; j < 4; j++)
        matrix[i][j] = 0.;

    for (i = 0; i < CS_DIM_3; i++)
      for (j = 0; j < CS_DIM_3; j++)
        matrix[i][j] = (double) mesh_builder->rotation_matrix[i][j];

    if (inv_pt == CS_TRUE) {

      for (i = 0; i < CS_DIM_3; i++)
        matrix[i][CS_DIM_3] = (double)
          - matrix[i][0]*mesh_builder->invariant_point[0]
          - matrix[i][1]*mesh_builder->invariant_point[1]
          - matrix[i][2]*mesh_builder->invariant_point[2]
          + mesh_builder->invariant_point[i];

      /* Clip "zero" values */

      for (i = 0; i < 3; i++)
        for (j = 0; j < 4; j++)
          if (FVM_ABS(matrix[i][j]) < 1.e-16)
            matrix[i][j] = 0.;

    }

    /* Now add periodicity */

    tr_id = fvm_periodicity_add_by_matrix(mesh->periodicity,
                                          mesh_builder->perio_num,
                                          FVM_PERIODICITY_ROTATION,
                                          matrix);

    bft_printf(" [ok]\n");

  } /* End if FVM_PERIODICITY_ROTATION */

}

/*============================================================================
 *  Public functions definition for API Fortran
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Receive messages from the pre-processor about the dimensions of mesh
 * parameters
 *
 * FORTRAN Interface:
 *
 * SUBROUTINE LEDEVI(NOMRUB, TYPENT, NBRENT, TABENT)
 * *****************
 *
 * INTEGER          NDIM        : <-- : Dimension de l'espace (3)
 * INTEGER          NCEL        : <-- : Nombre d'éléments actifs
 * INTEGER          NFAC        : <-- : Nombre de faces internes
 * INTEGER          NFABOR      : <-- : Nombre de faces de bord
 * INTEGER          NFML        : <-- : Nombre de familles des faces de bord
 * INTEGER          NPRFML      : <-- : Nombre de propriétés max par famille
 * INTEGER          NSOM        : <-- : Nombre de sommets (optionnel)
 * INTEGER          LNDFAC      : <-- : Longueur de SOMFAC (optionnel)
 * INTEGER          LNDFBR      : <-- : Longueur de SOMFBR (optionnel)
 * INTEGER          IPERIO      : <-- : Indicateur de périodicité
 * INTEGER          IPEROT      : <-- : Nombre de périodicités de rotation
 * INTEGER          ISTOP       : <-- : Indicateur d'arrêt demandé si > 0
 *----------------------------------------------------------------------------*/

void CS_PROCF(ledevi, LEDEVI)
(
 cs_int_t   *const ndim,    /* <-- dimension de l'espace                      */
 cs_int_t   *const ncel,    /* <-- nombre d'éléments actifs                   */
 cs_int_t   *const nfac,    /* <-- nombre de faces internes                   */
 cs_int_t   *const nfabor,  /* <-- nombre de faces de bord                    */
 cs_int_t   *const nfml,    /* <-- nombre de familles des faces de bord       */
 cs_int_t   *const nprfml,  /* <-- nombre de propriétés max par famille       */
 cs_int_t   *const nsom,    /* <-- nombre de sommets (optionnel)              */
 cs_int_t   *const lndfac,  /* <-- longueur de somfac (optionnel)             */
 cs_int_t   *const lndfbr,  /* <-- longueur de somfbr (optionnel)             */
 cs_int_t   *const iperio,  /* <-- indicateur de périodicité                  */
 cs_int_t   *const iperot   /* <-- nombre de périodicités de rotation         */
)
{
  cs_int_t  i;
  cs_pp_io_msg_header_t  header;

  cs_bool_t  dim_read = CS_FALSE;
  cs_bool_t  end_read = CS_FALSE;
  cs_pp_io_t  *comm = cs_glob_pp_io;
  cs_mesh_t  *mesh = cs_glob_mesh;
  cs_mesh_builder_t  *mesh_builder = cs_glob_mesh_builder;

  const char  *unexpected_msg = N_("Message de type <%s> sur <%s>\n"
                                  "inattendu ou de taille incorrecte");

  /* Initialize parameter values */

  *ndim = 3;
  *ncel = 0;
  *nfac = 0;
  *nfabor = 0;
  *nsom = 0;
  *lndfac = 0;
  *lndfbr = 0;
  *nfml = 0;
  *nprfml = 0;

  /* Loop onincomming messages */

  while (end_read == CS_FALSE) {

    /* Receive headers and clen header names */

    cs_pp_io_read_header(&header, comm);

    for (i = CS_PP_IO_NAME_LEN - 1; header.nom_rub[i] == ' '; i--);
    i++;
    if (i < CS_PP_IO_NAME_LEN) header.nom_rub[i] = '\0';

    /* Treatment according to the header name */

    if (strncmp(header.nom_rub, "EOF", CS_PP_IO_NAME_LEN)
        == 0) {
      cs_glob_pp_io = cs_pp_io_finalize(comm);
      comm = NULL;
    }

    if (strncmp(header.nom_rub, "start_block:dimensions",
                CS_PP_IO_NAME_LEN) == 0) {

      if (dim_read == CS_FALSE)
        dim_read = CS_TRUE;
      else
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.nom_rub, cs_pp_io_get_name(comm));

    }
    else if (strncmp(header.nom_rub, "end_block:dimensions",
                     CS_PP_IO_NAME_LEN) == 0) {

      if (dim_read == CS_TRUE) {
        dim_read = CS_FALSE;
        end_read = CS_TRUE;
      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.nom_rub, cs_pp_io_get_name(comm));

    }

    /* Receive dimensions from the pre-processor */

    else if (strncmp(header.nom_rub, "ndim",
                     CS_PP_IO_NAME_LEN) == 0) {

      if (dim_read != CS_TRUE || header.nbr_elt != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.nom_rub, cs_pp_io_get_name(comm));
      else
        cs_pp_io_read_body(&header, (void *) &(mesh->dim), comm);

    }
    else if (strncmp(header.nom_rub, "ncel",
                     CS_PP_IO_NAME_LEN) == 0) {

      if (dim_read != CS_TRUE || header.nbr_elt != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.nom_rub, cs_pp_io_get_name(comm));
      else {
        cs_pp_io_read_body(&header, (void *) &(mesh->n_cells), comm);
        mesh->n_g_cells = (fvm_gnum_t)mesh->n_cells;
      }

    }
    else if (strncmp(header.nom_rub, "nfac",
                     CS_PP_IO_NAME_LEN) == 0) {

      if (dim_read != CS_TRUE || header.nbr_elt != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.nom_rub, cs_pp_io_get_name(comm));
      else {
        cs_pp_io_read_body(&header, (void *) &(mesh->n_i_faces), comm);
        mesh->n_g_i_faces = (fvm_gnum_t)mesh->n_i_faces;
      }

    }
    else if (strncmp(header.nom_rub, "nfabor",
                     CS_PP_IO_NAME_LEN) == 0) {

      if (dim_read != CS_TRUE || header.nbr_elt != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.nom_rub, cs_pp_io_get_name(comm));
      else
        cs_pp_io_read_body(&header, (void *) &(mesh->n_b_faces), comm);

      mesh->n_g_b_faces = (fvm_gnum_t)mesh->n_b_faces;

    }
    else if (strncmp(header.nom_rub, "maillage:dim:familles:nbr",
                     CS_PP_IO_NAME_LEN) == 0) {
      if (dim_read != CS_TRUE || header.nbr_elt != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.nom_rub, cs_pp_io_get_name(comm));
      else
        cs_pp_io_read_body(&header, (void *) &(mesh->n_families), comm);

    }
    else if (strncmp(header.nom_rub, "maillage:dim:familles:prop:nbr",
                     CS_PP_IO_NAME_LEN) == 0) {

      if (dim_read != CS_TRUE || header.nbr_elt != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.nom_rub, cs_pp_io_get_name(comm));
      else
        cs_pp_io_read_body(&header,
                           (void *) &(mesh->n_max_family_items), comm);

    }
    else if (strncmp (header.nom_rub, "maillage:dim:groupes:nbr",
                      CS_PP_IO_NAME_LEN) == 0) {

      if (dim_read != CS_TRUE || header.nbr_elt != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.nom_rub, cs_pp_io_get_name(comm));
      else
        cs_pp_io_read_body(&header, (void *) &(mesh->n_groups), comm);

    }
    else if (strncmp(header.nom_rub, "maillage:dim:nsom",
                     CS_PP_IO_NAME_LEN) == 0) {

      if (dim_read != CS_TRUE || header.nbr_elt != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.nom_rub, cs_pp_io_get_name(comm));
      else {
        cs_pp_io_read_body(&header,
                           (void *) &(mesh->n_vertices), comm);
        mesh->n_g_vertices = mesh->n_vertices;
      }

    }
    else if (strncmp(header.nom_rub, "lndfac",
                     CS_PP_IO_NAME_LEN) == 0) {

      if (dim_read != CS_TRUE || header.nbr_elt != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.nom_rub, cs_pp_io_get_name(comm));
      else
        cs_pp_io_read_body(&header,
                           (void *) &(mesh->i_face_vtx_connect_size), comm);

    }
    else if (strncmp(header.nom_rub, "lndfbr",
                     CS_PP_IO_NAME_LEN) == 0) {

      if (dim_read != CS_TRUE || header.nbr_elt != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.nom_rub, cs_pp_io_get_name(comm));
      else
        cs_pp_io_read_body(&header,
                           (void *) &(mesh->b_face_vtx_connect_size), comm);

    }

    /* Additional messages for parallelism */

    else if (strncmp(header.nom_rub, "maillage:dim:nbrdom",
                     CS_PP_IO_NAME_LEN) == 0) {

      if (dim_read != CS_TRUE || header.nbr_elt != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.nom_rub, cs_pp_io_get_name(comm));
      else
        cs_pp_io_read_body(&header, (void *) &(mesh->n_domains),comm);

    }

    /* ----> Utilité de recevoir cette information ?
       = pour vérification ultérieure */

    else if (strncmp(header.nom_rub, "maillage:dim:inddom",
                     CS_PP_IO_NAME_LEN) == 0) {

      if (dim_read != CS_TRUE || header.nbr_elt != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.nom_rub, cs_pp_io_get_name(comm));
      else {
        cs_int_t  inddom;

        cs_pp_io_read_body(&header,
                           (void *) &inddom, comm);
        mesh->domain_num = inddom + 1;
      }

    }

    /* Additional messages for periodicity. Dimensions for periodic ghost
       cells have been received before. Here we allocate parameter list
       for periodicities and coupled face list for halo builder. */

    else if (strncmp(header.nom_rub, "maillage:dim:per:nbrtot",
                     CS_PP_IO_NAME_LEN) == 0) {

      if (dim_read != CS_TRUE || header.nbr_elt != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.nom_rub, cs_pp_io_get_name(comm));
      else {
        cs_pp_io_read_body(&header, (void *) &(mesh->n_init_perio), comm);

        assert(mesh->n_init_perio > 0);

        *iperio = 1;
        mesh->periodicity = fvm_periodicity_create(0.001);

        BFT_MALLOC(mesh_builder->per_face_idx,
                   mesh->n_init_perio + 1, cs_int_t);

      }

    }
    else if (strncmp(header.nom_rub, "maillage:dim:per:nbrrot",
                     CS_PP_IO_NAME_LEN) == 0) {

      if (dim_read != CS_TRUE || header.nbr_elt != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.nom_rub, cs_pp_io_get_name(comm));
      else
        cs_pp_io_read_body(&header, (void *) iperot, comm);

    }
    else if (strncmp(header.nom_rub, "maillage:dim:per:idxfac",
                     CS_PP_IO_NAME_LEN) == 0) {

      if (   dim_read != CS_TRUE
          || (cs_int_t)header.nbr_elt != mesh->n_init_perio + 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.nom_rub, cs_pp_io_get_name(comm));
      else
        cs_pp_io_read_body(&header,
                           (void *) &(mesh_builder->per_face_idx[0]),
                           comm);

      BFT_MALLOC(mesh_builder->per_face_lst,
                 2*mesh_builder->per_face_idx[mesh->n_init_perio], cs_int_t);

      if (cs_glob_base_nbr > 1)
        BFT_MALLOC(mesh_builder->per_rank_lst,
                   mesh_builder->per_face_idx[mesh->n_init_perio], cs_int_t);

    }
    else
      bft_error(__FILE__, __LINE__, 0,
                _(unexpected_msg), header.nom_rub, cs_pp_io_get_name(comm));

  } /* End of test on headers */

#if 0
  bft_printf("\nmesh->builder->per_face_idx:\n");
  for (i = 0; i < mesh->n_init_perio + 1; i++)
    bft_printf("per_face_idx[%d] = %d\n", i,
               mesh_builder->per_face_idx[i]);
  bft_printf_flush();
#endif

  /* Return values */

  *ndim = mesh->dim;
  *ncel = mesh->n_cells;
  *nfac = mesh->n_i_faces;
  *nfabor = mesh->n_b_faces;
  *nsom = mesh->n_vertices;
  *lndfac = mesh->i_face_vtx_connect_size;
  *lndfbr = mesh->b_face_vtx_connect_size;
  *nfml = mesh->n_families;
  *nprfml = mesh->n_max_family_items;

  /* Update data in cs_mesh_t structure in serial mode */

  if (cs_glob_base_nbr == 1) {

    mesh->n_cells_with_ghosts = mesh->n_cells;
    mesh->n_domains = 1;
    mesh->domain_num = 1;

  }

}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Receive data from the pre-processor and finalize communication with the
 * pre-processor
 *
 * mesh         <-- pointer to mesh structure
 *
 * returns:
 *----------------------------------------------------------------------------*/

void
cs_ecs_messages_read_data(cs_mesh_t      *const mesh)
{
  cs_int_t  i, perio_id, perio_num;
  cs_int_t  perio_type;
  cs_pp_io_msg_header_t  header;

  cs_int_t  n_elts = 0;
  cs_bool_t  end_read = CS_FALSE;
  cs_bool_t  data_read = CS_FALSE;
  cs_pp_io_t  *comm = cs_glob_pp_io;
  cs_mesh_builder_t  *mesh_builder = cs_glob_mesh_builder;

  const char  *unexpected_msg = N_("Message de type <%s> sur <%s>\n"
                                  "inattendu ou de taille incorrecte");

  /* Loop on incomming messages */

  while (end_read == CS_FALSE) {

    /* Receive header and clean header name */

    cs_pp_io_read_header(&header, comm);

    for (i = CS_PP_IO_NAME_LEN - 1 ; header.nom_rub[i] == ' ' ; i--);
    i++;
    if (i < CS_PP_IO_NAME_LEN) header.nom_rub[i] = '\0';

    /* Treatment according to the header name */

    if (strncmp(header.nom_rub, "EOF", CS_PP_IO_NAME_LEN)
        == 0) {
      cs_glob_pp_io = cs_pp_io_finalize (comm);
      comm = NULL;
    }

    if (strncmp(header.nom_rub, "start_block:data",
                CS_PP_IO_NAME_LEN) == 0) {

      if (data_read == CS_FALSE)
        data_read = CS_TRUE;
      else
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.nom_rub, cs_pp_io_get_name(comm));

    }
    else if (strncmp(header.nom_rub, "end_block:data",
                     CS_PP_IO_NAME_LEN) == 0) {

      if (data_read == CS_TRUE) {
        data_read = CS_FALSE;
        end_read = CS_TRUE;
      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.nom_rub, cs_pp_io_get_name(comm));

    }

    /* Receive data from the pre-processor */

    else if (strncmp(header.nom_rub, "ifacel",
                     CS_PP_IO_NAME_LEN) == 0) {

      n_elts = mesh->n_i_faces * 2;
      if (data_read != CS_TRUE || (cs_int_t)header.nbr_elt != n_elts)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.nom_rub, cs_pp_io_get_name(comm));
      else {
        BFT_MALLOC(mesh->i_face_cells, n_elts, cs_int_t);
        cs_pp_io_read_body(&header,
                           (void *) mesh->i_face_cells, comm);
      }

    }
    else if (strncmp(header.nom_rub, "ifabor",
                     CS_PP_IO_NAME_LEN) == 0) {

      if (data_read != CS_TRUE || (cs_int_t)header.nbr_elt != mesh->n_b_faces)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.nom_rub, cs_pp_io_get_name(comm));
      else {
        BFT_MALLOC(mesh->b_face_cells, mesh->n_b_faces, cs_int_t);
        cs_pp_io_read_body(&header,
                           (void *) mesh->b_face_cells, comm);
      }

    }
    else if (   strncmp(header.nom_rub, "maillage:data:fbr:famille",
                        CS_PP_IO_NAME_LEN) == 0
             || strncmp(header.nom_rub, "ifafbr",
                        CS_PP_IO_NAME_LEN) == 0) {

      if (data_read != CS_TRUE || (cs_int_t)header.nbr_elt != mesh->n_b_faces)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.nom_rub, cs_pp_io_get_name(comm));
      else {
        BFT_MALLOC(mesh->b_face_family, mesh->n_b_faces, cs_int_t);
        cs_pp_io_read_body(&header, (void *) mesh->b_face_family, comm);
      }

    }
    else if (strncmp(header.nom_rub, "maillage:data:cel:famille",
                     CS_PP_IO_NAME_LEN) == 0) {

      if (data_read != CS_TRUE || (cs_int_t)header.nbr_elt != mesh->n_cells)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.nom_rub, cs_pp_io_get_name(comm));
      else {
        /* Allocation for taking account of ghost cells */
        BFT_MALLOC(mesh->cell_family, mesh->n_cells, cs_int_t);
        cs_pp_io_read_body(&header, (void *) mesh->cell_family, comm);
      }

    }
    else if (   strncmp(header.nom_rub, "maillage:data:familles:propr",
                        CS_PP_IO_NAME_LEN) == 0
             || strncmp(header.nom_rub, "iprfml",
                        CS_PP_IO_NAME_LEN) == 0) {

      n_elts = mesh->n_families * mesh->n_max_family_items;
      if (data_read != CS_TRUE || (cs_int_t)header.nbr_elt != n_elts)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.nom_rub, cs_pp_io_get_name(comm));
      else {
        BFT_MALLOC(mesh->family_item, n_elts, cs_int_t);
        cs_pp_io_read_body(&header, (void *) mesh->family_item, comm);
      }

    }
    else if (strncmp(header.nom_rub, "maillage:data:groupes:pos",
                     CS_PP_IO_NAME_LEN) == 0) {

      if (data_read != CS_TRUE || (cs_int_t)header.nbr_elt != mesh->n_groups + 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.nom_rub, cs_pp_io_get_name(comm));
      else {
        BFT_MALLOC(mesh->group_idx, mesh->n_groups + 1, cs_int_t);
        cs_pp_io_read_body(&header, (void *) mesh->group_idx, comm);
      }

    }
    else if (strncmp(header.nom_rub, "maillage:data:groupes:noms",
                     CS_PP_IO_NAME_LEN) == 0) {

      if (data_read != CS_TRUE || mesh->group_idx == NULL ||
          (cs_int_t)header.nbr_elt != mesh->group_idx[mesh->n_groups] - 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.nom_rub, cs_pp_io_get_name(comm));
      else {
        BFT_MALLOC(mesh->group_lst, header.nbr_elt, char);
        cs_pp_io_read_body(&header, (void *) mesh->group_lst, comm);
      }

    }
    else if (   strncmp(header.nom_rub, "xyznod",
                        CS_PP_IO_NAME_LEN) == 0
             || strncmp(header.nom_rub, "maillage:data:som:xyz",
                        CS_PP_IO_NAME_LEN) == 0) {

      n_elts = mesh->dim * mesh->n_vertices;
      if (data_read != CS_TRUE || (cs_int_t)header.nbr_elt != n_elts)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.nom_rub, cs_pp_io_get_name(comm));
      else {
        BFT_MALLOC(mesh->vtx_coord, n_elts, cs_real_t);
        cs_pp_io_read_body(&header, (void *) mesh->vtx_coord, comm);
      }

    }
    else if (strncmp(header.nom_rub, "ipnfac",
                     CS_PP_IO_NAME_LEN) == 0) {

      n_elts = mesh->n_i_faces + 1;
      if (   data_read != CS_TRUE
          || mesh->n_vertices == 0 || (cs_int_t)header.nbr_elt != n_elts)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.nom_rub, cs_pp_io_get_name(comm));
      else {
        BFT_MALLOC(mesh->i_face_vtx_idx, n_elts, cs_int_t);
        cs_pp_io_read_body(&header,
                           (void *) mesh->i_face_vtx_idx, comm);
      }

    }
    else if (strncmp(header.nom_rub, "nodfac",
                     CS_PP_IO_NAME_LEN) == 0) {

      if (   data_read != CS_TRUE
          || mesh->n_vertices == 0
          || (cs_int_t)header.nbr_elt != mesh->i_face_vtx_connect_size)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.nom_rub, cs_pp_io_get_name(comm));
      else {
        BFT_MALLOC(mesh->i_face_vtx_lst, mesh->i_face_vtx_connect_size, cs_int_t);
        cs_pp_io_read_body(&header, (void *) mesh->i_face_vtx_lst, comm);
      }

    }
    else if (strncmp(header.nom_rub, "ipnfbr",
                     CS_PP_IO_NAME_LEN) == 0) {

      n_elts = mesh->n_b_faces + 1;
      if (   data_read != CS_TRUE
          || mesh->n_vertices == 0 || (cs_int_t)header.nbr_elt != n_elts)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.nom_rub, cs_pp_io_get_name(comm));
      else {
        BFT_MALLOC(mesh->b_face_vtx_idx, n_elts, cs_int_t);
        cs_pp_io_read_body(&header, (void *) mesh->b_face_vtx_idx, comm);
      }

    }
    else if (strncmp(header.nom_rub, "nodfbr",
                     CS_PP_IO_NAME_LEN) == 0) {

      if (   data_read != CS_TRUE
          || mesh->n_vertices == 0
          || (cs_int_t)header.nbr_elt != mesh->b_face_vtx_connect_size)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.nom_rub, cs_pp_io_get_name(comm));
      else {
        BFT_MALLOC(mesh->b_face_vtx_lst, mesh->b_face_vtx_connect_size, cs_int_t);
        cs_pp_io_read_body(&header, (void *) mesh->b_face_vtx_lst, comm);
      }

    }

    /* Additional buffers for parallelism */

    else if (strncmp(header.nom_rub, "maillage:data:cel:num",
                     CS_PP_IO_NAME_LEN) == 0) {

      if (data_read != CS_TRUE || (cs_int_t)header.nbr_elt != mesh->n_cells)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.nom_rub, cs_pp_io_get_name(comm));
      else {

        cs_int_t *recv_buf = NULL;

        BFT_MALLOC(mesh->global_cell_num, mesh->n_cells, fvm_gnum_t);
        BFT_MALLOC(recv_buf, mesh->n_cells, cs_int_t);
        cs_pp_io_read_body (&header, (void *) recv_buf, comm);
        for (i = 0; i < mesh->n_cells; i++)
          mesh->global_cell_num[i] = recv_buf[i];
        BFT_FREE(recv_buf);

      }
    }
    else if (strncmp(header.nom_rub, "maillage:data:fac:num",
                     CS_PP_IO_NAME_LEN) == 0) {

      if (data_read != CS_TRUE || (cs_int_t)header.nbr_elt != mesh->n_i_faces)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.nom_rub, cs_pp_io_get_name(comm));
      else {

        cs_int_t *recv_buf = NULL;

        BFT_MALLOC(mesh->global_i_face_num, mesh->n_i_faces, fvm_gnum_t);
        BFT_MALLOC(recv_buf, mesh->n_i_faces, cs_int_t);
        cs_pp_io_read_body(&header, (void *) recv_buf, comm);
        for (i = 0; i < mesh->n_i_faces; i++)
          mesh->global_i_face_num[i] = recv_buf[i];
        BFT_FREE(recv_buf);

      }
    }

    else if (strncmp(header.nom_rub, "maillage:data:fbr:num",
                     CS_PP_IO_NAME_LEN) == 0) {

      if (data_read != CS_TRUE || (cs_int_t)header.nbr_elt != mesh->n_b_faces)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.nom_rub, cs_pp_io_get_name(comm));
      else {

        cs_int_t  *recv_buf = NULL;

        BFT_MALLOC(mesh->global_b_face_num, mesh->n_b_faces, fvm_gnum_t);
        BFT_MALLOC(recv_buf, mesh->n_b_faces, cs_int_t);
        cs_pp_io_read_body(&header, (void *) recv_buf, comm);
        for (i = 0; i < mesh->n_b_faces; i++)
          mesh->global_b_face_num[i] = recv_buf[i];
        BFT_FREE(recv_buf);
      }

    }
    else if (strncmp(header.nom_rub, "maillage:data:som:num",
                     CS_PP_IO_NAME_LEN) == 0) {

      if (data_read != CS_TRUE || (cs_int_t)header.nbr_elt != mesh->n_vertices)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.nom_rub, cs_pp_io_get_name(comm));
      else {

        cs_int_t  *recv_buf = NULL;

        BFT_MALLOC(mesh->global_vtx_num, mesh->n_vertices, fvm_gnum_t);
        BFT_MALLOC(recv_buf, mesh->n_vertices, cs_int_t);
        cs_pp_io_read_body(&header, (void *) recv_buf, comm);
        for (i = 0; i < mesh->n_vertices; i++)
          mesh->global_vtx_num[i] = recv_buf[i];
        BFT_FREE(recv_buf);
      }

    }

    /* Additional buffers for periodicity */

    else if (strncmp(header.nom_rub, "maillage:data:per:type:",
                     strlen("maillage:data:per:type:")) == 0) {

      if (data_read != CS_TRUE || header.nbr_elt != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.nom_rub, cs_pp_io_get_name(comm));
      else {

        perio_num = atoi(header.nom_rub + strlen("maillage:data:per:type:"));
        mesh_builder->perio_num = perio_num;
        cs_pp_io_read_body(&header, (void *) &(perio_type), comm);
        mesh_builder->perio_type = (fvm_periodicity_type_t) perio_type;

      }

    }
    else if (strncmp(header.nom_rub, "maillage:data:per:trans:",
                     strlen("maillage:data:per:trans:")) == 0) {

      if (data_read != CS_TRUE || (cs_int_t)header.nbr_elt != mesh->dim)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.nom_rub, cs_pp_io_get_name(comm));
      else {

        perio_num = atoi(header.nom_rub + strlen("maillage:data:per:trans:"));
        assert(mesh_builder->perio_num == perio_num);
        cs_pp_io_read_body(&header, (void *) mesh_builder->translation, comm);

      }

    }
    else if (strncmp(header.nom_rub, "maillage:data:per:ptinv:",
                     strlen("maillage:data:per:ptinv:")) == 0) {

      if (data_read != CS_TRUE || (cs_int_t)header.nbr_elt != mesh->dim)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.nom_rub, cs_pp_io_get_name(comm));
      else {

        perio_num = atoi(header.nom_rub + strlen("maillage:data:per:ptinv:"));
        assert(mesh_builder->perio_num == perio_num);
        cs_pp_io_read_body(&header,
                           (void *) mesh_builder->invariant_point,
                           comm);
      }

    }
    else if (strncmp(header.nom_rub, "maillage:data:per:matrice:",
                     strlen("maillage:data:per:matrice:")) == 0) {

      n_elts = mesh->dim * mesh->dim;
      if (data_read != CS_TRUE || (cs_int_t)header.nbr_elt != n_elts)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.nom_rub, cs_pp_io_get_name(comm));
      else {

        perio_num = atoi(header.nom_rub + strlen("maillage:data:per:matrice:"));
        assert(mesh_builder->perio_num == perio_num);
        cs_pp_io_read_body(&header,
                           (void *) mesh_builder->rotation_matrix,
                           comm);

        /* Add a periodicity to mesh->periodicities */

        _add_periodicity();

      }

    }
    else if (strncmp(header.nom_rub, "maillage:data:per:fac:",
                     strlen("maillage:data:per:fac:")) == 0) {

      perio_id = atoi(header.nom_rub
                     + strlen("maillage:data:per:fac:")) - 1;
      n_elts = 2 * (  mesh_builder->per_face_idx[perio_id + 1]
                    - mesh_builder->per_face_idx[perio_id] );

      if (data_read != CS_TRUE || (cs_int_t)header.nbr_elt != n_elts)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.nom_rub, cs_pp_io_get_name(comm));
      else {

        cs_int_t  shift = 2 * mesh_builder->per_face_idx[perio_id];
        cs_int_t  *per_face_lst = mesh_builder->per_face_lst + shift;

        cs_pp_io_read_body(&header,
                           (void *)per_face_lst,
                           comm);

      }

    }
    else if (strncmp(header.nom_rub, "maillage:data:per:dom:",
                     strlen("maillage:data:per:dom:")) == 0) {

      perio_id = atoi(header.nom_rub
                     + strlen("maillage:data:per:dom:")) - 1;
      n_elts =  mesh_builder->per_face_idx[perio_id + 1]
              - mesh_builder->per_face_idx[perio_id];

      if (data_read != CS_TRUE || (cs_int_t)header.nbr_elt != n_elts)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.nom_rub, cs_pp_io_get_name(comm));
      else {

        cs_int_t  shift = mesh_builder->per_face_idx[perio_id];
        cs_int_t  *dom_fac_perio = mesh_builder->per_rank_lst + shift;

        cs_pp_io_read_body(&header,
                           (void *)dom_fac_perio,
                           comm);

      }

    }

  } /* End of loop on messages */

  /* Finalize communications with the pre-processor */
  /*------------------------------------------------*/

  if (cs_glob_pp_io != NULL) {
    cs_pp_io_finalize(cs_glob_pp_io);
    cs_glob_pp_io = NULL;
  }

  /* Define a new vertex numbering in parallel mode because the global
     numbering from the pre-processor doesn't insure an increasing
     numbering */

  if (mesh->global_vtx_num != NULL)
    cs_mesh_order_vertices(mesh);

#if 0 /* For debugging purpose */
  bft_printf("\nmesh_builder->per_face_lst && mesh_builder->per_rank_lst:\n");
  for (perio_id = 0; perio_id < mesh->n_init_perio; perio_id++) {
    for (i = mesh_builder->per_face_idx[perio_id];
         i < mesh_builder->per_face_idx[perio_id+1]; i++) {
      if (cs_glob_base_nbr > 1)
        bft_printf("\n\tICOUPL: %d, per_rank_lst: %d\n",
                   i, mesh_builder->per_rank_lst[i]);
      bft_printf("\tper_face_lst[%d] = %d\n",
                 2*i, mesh_builder->per_face_lst[2*i]);
      bft_printf("\tper_face_lst[%d] = %d\n",
                 2*i+1, mesh_builder->per_face_lst[2*i+1]);
    }
  }
  bft_printf_flush();
#endif

}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

