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

#ifndef __CS_HALO_H__
#define __CS_HALO_H__

/*============================================================================
 * Structure and function headers handling with ghost cells
 *============================================================================*/

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <fvm_interface.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/* Halo type */

typedef enum {

  CS_HALO_STANDARD,
  CS_HALO_EXTENDED,
  CS_HALO_N_TYPES

} cs_halo_type_t;

/*============================================================================
 * Type definition
 *============================================================================*/

/* Structure for halo management */
/* ----------------------------- */

typedef struct {

  cs_int_t  n_c_domains;     /* Number of communicating domains. */
  cs_int_t  *c_domain_rank;  /* List of communicating ranks */

  /* send_halo features : send to distant ranks */

  cs_int_t  n_send_elts[2];    /* Numer of ghost elements in send_halo
                                n_elts[0] = standard elements
                                n_elts[1] = extended + standard elements */

  cs_int_t  *send_list;        /* List of local numbers of elements in send_halo */

  cs_int_t  *send_index;       /* Index on in_elements.
                                  Size = 2*n_c_domains. For each rank, we
                                  have an index for standard halo and one
                                  for extended halo. */

  cs_int_t  *send_perio_lst ;  /* For each transformation and for each type of halo
                                  on each communicating rank, we store 2 data:
                                   - start index,
                                   - number of elements. */

  /* halo features : receive from distant ranks */

  cs_int_t  n_elts[2];      /* Numer of ghost elements in halo
                                 n_elts[0] = standard elements
                                 n_elts[1] = extended + standard elements */

  cs_int_t  *list;          /* List of local numbers of elements in halo */

  cs_int_t  *index;         /* Index on in_elements.
                               Size = 2*n_c_domains. For each rank, we
                               have an index for standard halo and one
                               for extended halo. */

  cs_int_t  *perio_lst;     /* For each transformation and for each type of halo
                               on each communicating rank, we store 2 data:
                                 - start index,
                                 - number of elements. */

  /* Variables used during the synchronization process */

  cs_real_t  *tmp_buffer;   /* Buffer used to de-interlace variable
                               in case of strided variable to sync. */

#if defined(_CS_HAVE_MPI)
  MPI_Request   *mpi_request;   /* MPI Request array */
  MPI_Status    *mpi_status;    /* MPI Status array */

  cs_real_t  *comm_buffer;      /* Buffer for the communication purpose.
                                   Buffer size is equal to the maximum
                                   number of ghost cells between send_halo and
                                   halo. */
#endif

  /* Organisation of perio_lst:

         -------------------------------------------------
    T1:  |   |   |   |   |   |   |   |   |   |   |   |   |
         -------------------------------------------------
          idx  n  idx  n  idx  n  idx  n  idx  n  idx  n
          ______  ______  ______  ______  ______  ______
           std     ext     std     ext     std     ext
           ___________     ___________     ___________
             rank 0          rank 1          rank 2

         -------------------------------------------------
    T2:  |   |   |   |   |   |   |   |   |   |   |   |   |
         -------------------------------------------------
          idx  n  idx  n  idx  n  idx  n  idx  n  idx  n
          ______  ______  ______  ______  ______  ______
           std     ext     std     ext     std     ext
           ___________     ___________     ___________
             rank 0          rank 1          rank 2

         -------------------------------------------------
    T3:  |   |   |   |   |   |   |   |   |   |   |   |   |
         -------------------------------------------------
          idx  n  idx  n  idx  n  idx  n  idx  n  idx  n
          ______  ______  ______  ______  ______  ______
           std     ext     std     ext     std     ext
           ___________     ___________     ___________
             rank 0          rank 1          rank 2

  etc...

  */

} cs_halo_t;

/*=============================================================================
 * Global static variables
 *============================================================================*/

/*============================================================================
 *  Public function header for Fortran API
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create a halo structure.
 *
 * parameters:
 *   ifs  -->  pointer to a fvm_interface_set structure
 *
 * returns:
 *   pointer to created cs_mesh_halo_t structure
 *---------------------------------------------------------------------------*/

cs_halo_t *
cs_halo_create(fvm_interface_set_t  *ifs);

/*----------------------------------------------------------------------------
 * Destroy a halo structure.
 *
 * parameters:
 *   this_halo  -->  pointer to cs_mesh_halo structure to destroy
 *
 * returns:
 *   pointer to deleted halo structure (NULL)
 *---------------------------------------------------------------------------*/

cs_halo_t *
cs_halo_destroy(cs_halo_t  *this_halo);

/*----------------------------------------------------------------------------
 * Dump a cs_mesh_halo_t structure.
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
             cs_halo_t  *halo);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CS_HALO_H__ */
