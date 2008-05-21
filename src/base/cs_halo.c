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
 * Functions dealing with ghost cells
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <string.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_error.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <fvm_interface.h>
#include <fvm_periodicity.h>
#include <fvm_order.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_halo.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create a halo structure.
 *
 * parameters:
 *   ifs  -->  pointer to a fvm_interface_set structure
 *
 * returns:
 *  pointer to created cs_mesh_halo_t structure
 *---------------------------------------------------------------------------*/

cs_halo_t *
cs_halo_create(fvm_interface_set_t  *ifs)
{
  cs_int_t  i, tmp_id, perio_lst_size;

  cs_int_t  n_transforms = 0;
  cs_int_t  loc_id = -1;

  cs_halo_t  *halo = NULL;

  const fvm_interface_t  *interface = NULL;
  const fvm_periodicity_t  *periodicity = fvm_interface_set_periodicity(ifs);

  BFT_MALLOC(halo, 1, cs_halo_t);

  halo->n_c_domains = fvm_interface_set_size(ifs);
  BFT_MALLOC(halo->c_domain_rank, halo->n_c_domains, cs_int_t);

  /* Check if cs_glob_base_rang belongs to interface set in order to
     order ranks with local rank at first place */

  for (i = 0; i < halo->n_c_domains; i++) {

    interface = fvm_interface_set_get(ifs, i);
    halo->c_domain_rank[i] = fvm_interface_rank(interface);

    if (cs_glob_base_rang == fvm_interface_rank(interface))
      loc_id = i;

  } /* End of loop on ranks */

  if (loc_id > 0) {

    tmp_id = halo->c_domain_rank[loc_id];
    halo->c_domain_rank[loc_id] = halo->c_domain_rank[0];
    halo->c_domain_rank[0] = tmp_id;

  }

  /* Order ranks */

  if (   halo->n_c_domains > 2
      && fvm_order_local_test(&(halo->c_domain_rank[1]),
                              NULL,
                              halo->n_c_domains-1) == 0) {

    fvm_lnum_t  *order = NULL;
    fvm_gnum_t  *buffer = NULL;

    assert(sizeof(fvm_lnum_t) == sizeof(cs_int_t));

    BFT_MALLOC(order, halo->n_c_domains - 1, fvm_lnum_t);
    BFT_MALLOC(buffer, halo->n_c_domains - 1, fvm_gnum_t);

    for (i = 1; i < halo->n_c_domains; i++)
      buffer[i-1] = (fvm_gnum_t)halo->c_domain_rank[i];

    fvm_order_local_allocated(NULL,
                              buffer,
                              order,
                              halo->n_c_domains - 1);

    for (i = 0; i < halo->n_c_domains - 1; i++)
      halo->c_domain_rank[i+1] = (cs_int_t)buffer[order[i]];

    BFT_FREE(buffer);
    BFT_FREE(order);

  } /* End of ordering ranks */

  n_transforms = fvm_periodicity_get_n_transforms(periodicity);

  /* We need 2 data per transformation and there are n_transforms
     transformations. For each rank, we need data for standard and
     extended halo. */

  perio_lst_size = 2*n_transforms * 2*halo->n_c_domains;

  BFT_MALLOC(halo->send_perio_lst, perio_lst_size, cs_int_t);
  BFT_MALLOC(halo->perio_lst, perio_lst_size, cs_int_t);

  for (i = 0; i < perio_lst_size; i++) {
    halo->send_perio_lst[i] = 0;
    halo->perio_lst[i] = 0;
  }

  BFT_MALLOC(halo->send_index, 2*halo->n_c_domains + 1, cs_int_t);
  BFT_MALLOC(halo->index, 2*halo->n_c_domains + 1, cs_int_t);

  for (i = 0; i < 2*halo->n_c_domains + 1; i++) {
    halo->send_index[i] = 0;
    halo->index[i] = 0;
  }

  halo->send_list = NULL;
  halo->list = NULL;
  halo->tmp_buffer = NULL;

#if defined(_CS_HAVE_MPI)
  BFT_MALLOC(halo->mpi_request, 2*halo->n_c_domains, MPI_Request);
  BFT_MALLOC(halo->mpi_status, 2*halo->n_c_domains, MPI_Status);

  halo->comm_buffer = NULL;
#endif

  return halo;
}

/*----------------------------------------------------------------------------
 * Destroy a halo structure
 *
 * parameters:
 *   halo  -->  pointer to cs_mesh_halo_t structure to destroy
 *
 * Returns:
 *  pointer to deleted halo structure (NULL)
 *---------------------------------------------------------------------------*/

cs_halo_t *
cs_halo_destroy(cs_halo_t  *halo)
{
  if (halo == NULL)
    return NULL;

  halo->n_c_domains = 0;
  BFT_FREE(halo->c_domain_rank);

  BFT_FREE(halo->send_perio_lst);
  BFT_FREE(halo->send_index);
  BFT_FREE(halo->perio_lst);
  BFT_FREE(halo->index);

  if (halo->send_list != NULL)
    BFT_FREE(halo->send_list);

  if (halo->list != NULL)
    BFT_FREE(halo->list);

  if (halo->tmp_buffer != NULL)
    BFT_FREE(halo->tmp_buffer);

#if defined(_CS_HAVE_MPI)
  BFT_FREE(halo->mpi_request);
  BFT_FREE(halo->mpi_status);
  BFT_FREE(halo->comm_buffer);
#endif

  BFT_FREE(halo);

  return NULL;
}

/*----------------------------------------------------------------------------
 * Dump a cs_halo_t structure.
 *
 * parameters:
 *   n_cells        -->  number of cells
 *   n_init_perio   -->  initial number of periodicity
 *   n_transforms   -->  number of transformations
 *   print_level    -->  0 only dimensions and indexes are printed, else (1)
 *                       everything is printed
 *   halo           --> pointer to cs_halo_t struture
 *---------------------------------------------------------------------------*/

void
cs_halo_dump(cs_int_t    n_cells,
             cs_int_t    n_init_perio,
             cs_int_t    n_transforms,
             cs_int_t    print_level,
             cs_halo_t  *halo)
{
  cs_int_t  i, j, halo_id;

  if (halo == NULL) {
    bft_printf(_("\n\n  halo: nil\n"));
    return;
  }

  bft_printf(_("\n  halo        : %p\n"
               "  n_init_perio  : %d\n"
               "  n_transforms  : %d\n"
               "  n_c_domains   : %d\n"),
             halo, n_init_perio, n_transforms,
             halo->n_c_domains);

  bft_printf("\nRanks on halo frontier :\n");
  for (i = 0; i < halo->n_c_domains; i++)
    bft_printf("%5d", halo->c_domain_rank[i]);

  for (halo_id = 0; halo_id < 2; halo_id++) {

    cs_int_t  n_elts[2];
    cs_int_t  *index = NULL, *list = NULL, *perio_lst = NULL;

    bft_printf("\n    ---------\n");

    if (halo_id == 0) {

      bft_printf("    send_halo :\n");
      n_elts[0] = halo->n_send_elts[0];
      n_elts[1] = halo->n_send_elts[1];
      index = halo->send_index;
      list = halo->send_list;
      perio_lst = halo->send_perio_lst;

    }
    else if (halo_id == 1) {

      bft_printf("    halo :\n");
      n_elts[0] = halo->n_elts[0];
      n_elts[1] = halo->n_elts[1];
      index = halo->index;
      list = halo->list;
      perio_lst = halo->perio_lst;

    }

    bft_printf("    ---------\n\n");
    bft_printf(_("  n_ghost_cells       : %d\n"
                 "  n_std_ghost_cells   : %d\n"), n_elts[1], n_elts[0]);

    if (index == NULL)
      return;

    if (n_init_perio > 0) {

      const cs_int_t  stride = 4*halo->n_c_domains;

      for (i = 0; i < n_transforms; i++) {

        bft_printf("\nTransformation n°: %d\n", i+1);

        for (j = 0; j < halo->n_c_domains; j++) {

          bft_printf("    rank %3d <STD> %5d %5d <EXT> %5d %5d\n",
                     halo->c_domain_rank[j],
                     perio_lst[i*stride + 4*j],
                     perio_lst[i*stride + 4*j+1],
                     perio_lst[i*stride + 4*j+2],
                     perio_lst[i*stride + 4*j+3]);
        }

      } /* End of loop on perio */

    } /* End if n_perio > 0 */

    for (i = 0; i < halo->n_c_domains; i++) {

      bft_printf(_("\n  rank      %d:\n"), halo->c_domain_rank[i]);

      if (index[2*i+1] - index[2*i] > 0) {

        bft_printf(_("\n  Standard halo\n"));
        bft_printf(_("  idx start %d:          idx end   %d:\n"),
                   index[2*i], index[2*i+1]);

        if (print_level == 1) {
          bft_printf(_("\n            id      cell number\n"));
          for (j = index[2*i]; j < index[2*i+1]; j++)
            bft_printf(_("    %10d %10d %10d\n"),
                       j, list[j]+1, n_cells+j+1);
        }

      } /* there are elements on standard neighborhood */

      if (index[2*i+2] - index[2*i+1] > 0) {

        bft_printf(_("\n  Extended halo\n"));
        bft_printf(_("  idx start %d:          idx end   %d:\n"),
                   index[2*i+1], index[2*i+2]);

        if (print_level == 1) {
          bft_printf(_("\n            id      cell number\n"));
          for (j = index[2*i+1]; j < index[2*i+2]; j++)
            bft_printf(_("    %10d %10d %10d\n"),
                       j, list[j]+1, n_cells+j+1);
        }

      } /* If there are elements on extended neighborhood */

    } /* End of loop on involved ranks */

  } /* End of loop on halos (send_halo/halo) */

  bft_printf("\n\n");
  bft_printf_flush();
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
