#ifndef __CS_HALO_H__
#define __CS_HALO_H__

/*============================================================================
 * Structure and function headers handling with ghost cells
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_interface.h"

#include "fvm_periodicity.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Type definitions
 *============================================================================*/

/* Halo type */

typedef enum {

  CS_HALO_STANDARD,
  CS_HALO_EXTENDED,
  CS_HALO_N_TYPES

} cs_halo_type_t;

/* Options for cs_halo_sync_component(). */

typedef enum {

  CS_HALO_ROTATION_COPY,     /* Copy halo */
  CS_HALO_ROTATION_ZERO,     /* Set rotation halo values to zero */
  CS_HALO_ROTATION_IGNORE    /* Do not modify rotation halo values */

} cs_halo_rotation_t ;

/* Structure for halo management */
/* ----------------------------- */

typedef struct {

  int       n_c_domains;     /* Number of communicating domains. */
  int       n_transforms;    /* Number of periodic transformations */

  int       *c_domain_rank;  /* List of communicating ranks */

  const fvm_periodicity_t * periodicity; /* Pointer to periodicity
                                            structure describing transforms */

  int       n_rotations;     /* Number of periodic transformations
                                involving rotations */

  cs_lnum_t  n_local_elts;   /* Number of local elements */

  /* send_halo features : send to distant ranks */

  cs_lnum_t  n_send_elts[2];   /* Numer of ghost elements in send_list
                                n_elts[0] = standard elements
                                n_elts[1] = extended + standard elements */

  cs_lnum_t  *send_list;       /* List of local elements in distant halos
                                  (0 to n-1 numbering) */

  cs_lnum_t  *send_index;      /* Index on send_list
                                  Size = 2*n_c_domains + 1. For each rank, we
                                  have an index for standard halo and one
                                  for extended halo. */

  cs_lnum_t  *send_perio_lst ; /* For each transformation and for each type of
                                  halo on each communicating rank, we store
                                  2 values:
                                   - start index,
                                   - number of elements. */

  /* halo features : receive from distant ranks */

  cs_lnum_t  n_elts[2];       /* Numer of ghost elements in halo
                                 n_elts[0] = standard elements
                                 n_elts[1] = extended + standard elements */

  cs_lnum_t  *index;        /* Index on halo sections;
                               Size = 2*n_c_domains. For each rank, we
                               have an index for the standard halo and one
                               for the extended halo. */

  cs_lnum_t  *perio_lst;    /* For each transformation and for each type of halo
                               on each communicating rank, we store 2 values:
                                 - start index,
                                 - number of elements. */

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
 *   ifs  <--  pointer to a fvm_interface_set structure
 *
 * returns:
 *   pointer to created cs_halo_t structure
 *---------------------------------------------------------------------------*/

cs_halo_t *
cs_halo_create(cs_interface_set_t  *ifs);

/*----------------------------------------------------------------------------
 * Create a halo structure, using a reference halo
 *
 * parameters:
 *   ref  <--  pointer to reference halo
 *
 * returns:
 *   pointer to created cs_halo_t structure
 *---------------------------------------------------------------------------*/

cs_halo_t *
cs_halo_create_from_ref(const cs_halo_t  *ref);

/*----------------------------------------------------------------------------
 * Destroy a halo structure.
 *
 * parameters:
 *   this_halo  <--  pointer to cs_halo structure to destroy
 *
 * returns:
 *   pointer to deleted halo structure (NULL)
 *---------------------------------------------------------------------------*/

cs_halo_t *
cs_halo_destroy(cs_halo_t  *this_halo);

/*----------------------------------------------------------------------------
 * Update global buffer sizes so as to be usable with a given halo.
 *
 * Calls to halo synchronizations with variable strides up to 3 are
 * expected. For strides greater than 3, the halo will be resized if
 * necessary directly by the synchronization function.
 *
 * This function should be called at the end of any halo creation,
 * so that buffer sizes are increased if necessary.
 *
 * parameters:
 *   halo <-- pointer to cs_halo_t structure.
 *---------------------------------------------------------------------------*/

void
cs_halo_update_buffers(const cs_halo_t *halo);

/*----------------------------------------------------------------------------
 * Free global halo backup buffer.
 *---------------------------------------------------------------------------*/

void
cs_halo_free_buffer(void);

/*----------------------------------------------------------------------------
 * Apply local cells renumbering to a halo
 *
 * parameters:
 *   halo        <-- pointer to halo structure
 *   new_cell_id <-- array indicating old -> new cell id (0 to n-1)
 *---------------------------------------------------------------------------*/

void
cs_halo_renumber_cells(cs_halo_t        *halo,
                       const cs_lnum_t   new_cell_id[]);

/*----------------------------------------------------------------------------
 * Apply ghost cells renumbering to a halo
 *
 * parameters:
 *   halo        <-- pointer to halo structure
 *   old_cell_id <-- array indicating new -> old cell id (0 to n-1)
 *---------------------------------------------------------------------------*/

void
cs_halo_renumber_ghost_cells(cs_halo_t        *halo,
                             const cs_lnum_t   old_cell_id[]);

/*----------------------------------------------------------------------------
 * Update array of any type of halo values in case of parallelism or
 * periodicity.
 *
 * Data is untyped; only its size is given, so this function may also
 * be used to synchronize interleaved multidimendsional data, using
 * size = element_size*dim (assuming a homogeneous environment, at least
 * as far as data encoding goes).
 *
 * This function aims at copying main values from local elements
 * (id between 1 and n_local_elements) to ghost elements on distant ranks
 * (id between n_local_elements + 1 to n_local_elements_with_halo).
 *
 * parameters:
 *   halo      <-- pointer to halo structure
 *   sync_mode <-- synchronization mode (standard or extended)
 *   num       <-> pointer to local number value array
 *----------------------------------------------------------------------------*/

void
cs_halo_sync_untyped(const cs_halo_t  *halo,
                     cs_halo_type_t    sync_mode,
                     size_t            size,
                     void             *val);

/*----------------------------------------------------------------------------
 * Update array of integer halo values in case of parallelism or periodicity.
 *
 * This function aims at copying main values from local elements
 * (id between 1 and n_local_elements) to ghost elements on distant ranks
 * (id between n_local_elements + 1 to n_local_elements_with_halo).
 *
 * parameters:
 *   halo      <-- pointer to halo structure
 *   sync_mode <-- synchronization mode (standard or extended)
 *   num       <-> pointer to local number value array
 *----------------------------------------------------------------------------*/

void
cs_halo_sync_num(const cs_halo_t  *halo,
                 cs_halo_type_t    sync_mode,
                 cs_lnum_t         num[]);

/*----------------------------------------------------------------------------
 * Update array of variable (floating-point) halo values in case of
 * parallelism or periodicity.
 *
 * This function aims at copying main values from local elements
 * (id between 1 and n_local_elements) to ghost elements on distant ranks
 * (id between n_local_elements + 1 to n_local_elements_with_halo).
 *
 * parameters:
 *   halo      <-- pointer to halo structure
 *   sync_mode <-- synchronization mode (standard or extended)
 *   var       <-> pointer to variable value array
 *----------------------------------------------------------------------------*/

void
cs_halo_sync_var(const cs_halo_t  *halo,
                 cs_halo_type_t    sync_mode,
                 cs_real_t         var[]);

/*----------------------------------------------------------------------------
 * Update array of strided variable (floating-point) halo values in case
 * of parallelism or periodicity.
 *
 * This function aims at copying main values from local elements
 * (id between 1 and n_local_elements) to ghost elements on distant ranks
 * (id between n_local_elements + 1 to n_local_elements_with_halo).
 *
 * parameters:
 *   halo      <-- pointer to halo structure
 *   sync_mode <-- synchronization mode (standard or extended)
 *   var       <-> pointer to variable value array
 *   stride    <-- number of (interlaced) values by entity
 *----------------------------------------------------------------------------*/

void
cs_halo_sync_var_strided(const cs_halo_t  *halo,
                         cs_halo_type_t    sync_mode,
                         cs_real_t         var[],
                         int               stride);

/*----------------------------------------------------------------------------
 * Update array of vector variable component (floating-point) halo values
 * in case of parallelism or periodicity.
 *
 * This function aims at copying main values from local elements
 * (id between 1 and n_local_elements) to ghost elements on distant ranks
 * (id between n_local_elements + 1 to n_local_elements_with_halo).
 *
 * If rotation_op is equal to CS_HALO_ROTATION_IGNORE, halo values
 * corresponding to periodicity with rotation are left unchanged from their
 * previous values.
 *
 * If rotation_op is equal to CS_HALO_ROTATION_ZERO, halo values
 * corresponding to periodicity with rotation are set to 0.
 *
 * If rotation_op is equal to CS_HALO_ROTATION_COPY, halo values
 * corresponding to periodicity with rotation are exchanged normally, so
 * the behavior is the same as that of cs_halo_sync_var().
 *
 * parameters:
 *   halo        <-- pointer to halo structure
 *   sync_mode   <-- synchronization mode (standard or extended)
 *   rotation_op <-- rotation operation
 *   var         <-> pointer to variable value array
 *----------------------------------------------------------------------------*/

void
cs_halo_sync_component(const cs_halo_t    *halo,
                       cs_halo_type_t      sync_mode,
                       cs_halo_rotation_t  rotation_op,
                       cs_real_t           var[]);

/*----------------------------------------------------------------------------
 * Update array of strided vector variable components (floating-point)
 * halo values in case of parallelism or periodicity.
 *
 * This function aims at copying main values from local elements
 * (id between 1 and n_local_elements) to ghost elements on distant ranks
 * (id between n_local_elements + 1 to n_local_elements_with_halo).
 *
 * If rotation_op is equal to CS_HALO_ROTATION_IGNORE, halo values
 * corresponding to periodicity with rotation are left unchanged from their
 * previous values.
 *
 * If rotation_op is equal to CS_HALO_ROTATION_ZERO, halo values
 * corresponding to periodicity with rotation are set to 0.
 *
 * If rotation_op is equal to CS_HALO_ROTATION_COPY, halo values
 * corresponding to periodicity with rotation are exchanged normally, so
 * the behavior is the same as that of cs_halo_sync_var_strided().
 *
 * parameters:
 *   halo        <-- pointer to halo structure
 *   sync_mode   <-- synchronization mode (standard or extended)
 *   rotation_op <-- rotation operation
 *   var         <-> pointer to variable value array
 *   stride      <-- number of (interlaced) values by entity
 *----------------------------------------------------------------------------*/

void
cs_halo_sync_components_strided(const cs_halo_t    *halo,
                                cs_halo_type_t      sync_mode,
                                cs_halo_rotation_t  rotation_op,
                                cs_real_t           var[],
                                int                 stride);

/*----------------------------------------------------------------------------
 * Return MPI_Barrier usage flag.
 *
 * returns:
 *   true if MPI barriers are used after posting receives and before posting
 *   sends, false otherwise
 *---------------------------------------------------------------------------*/

bool
cs_halo_get_use_barrier(void);

/*----------------------------------------------------------------------------
 * Set MPI_Barrier usage flag.
 *
 * parameters:
 *   use_barrier <-- true if MPI barriers should be used after posting
 *                   receives and before posting sends, false otherwise.
 *---------------------------------------------------------------------------*/

void
cs_halo_set_use_barrier(bool use_barrier);

/*----------------------------------------------------------------------------
 * Dump a cs_halo_t structure.
 *
 * parameters:
 *   halo           <--  pointer to cs_halo_t struture
 *   print_level    <--  0 only dimensions and indexes are printed, else (1)
 *                       everything is printed
 *---------------------------------------------------------------------------*/

void
cs_halo_dump(const cs_halo_t  *halo,
             int               print_level);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_HALO_H__ */
