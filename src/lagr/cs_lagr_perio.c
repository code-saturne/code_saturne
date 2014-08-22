/*============================================================================
 * Management of the periodicity for particles
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2014 EDF S.A.

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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_printf.h"

#include "fvm_periodicity.h"

#include "cs_base.h"
#include "cs_halo.h"
#include "cs_mesh.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_perio.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 *  Public functions definition for API Fortran
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Build buffers to keep the link between a given halo cell and :
 *  - the related real cell
 *  - the transformation id
 *
 * Fortran Interface :
 *
 * SUBROUTINE PERLOC
 * *****************
 *
 * INTEGER ICELCR(NCELET-NCEL) : <-  : related real cell buffer
 * INTEGER IPERCR(NCELET-NCEL) : <-  : transformation id buffer
 *
 * Returns:
 *----------------------------------------------------------------------------*/

void
CS_PROCF (perloc, PERLOC)(cs_int_t   *icelcr,
                          cs_int_t   *ipercr)
{
  cs_int_t  i, rank_id, shift, t_id;
  cs_int_t  start_std, end_std, length, start_ext, end_ext;

  cs_mesh_t  *mesh = cs_glob_mesh;
  cs_halo_t  *halo = mesh->halo;

  const cs_int_t  n_transforms = mesh->n_transforms;
  const cs_int_t  local_rank = (cs_glob_rank_id == -1) ? 0:cs_glob_rank_id;

  if (mesh->halo_type == CS_HALO_N_TYPES)
    return;

  assert(halo != NULL);

  for (t_id = 0; t_id < n_transforms; t_id++) {

    shift = 4 * halo->n_c_domains * t_id;

    for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

      if (   mesh->n_domains == 1
          || halo->c_domain_rank[rank_id] == local_rank) {

        start_std = halo->perio_lst[shift + 4*rank_id];
        length = halo->perio_lst[shift + 4*rank_id + 1];
        end_std = start_std + length;

        for (i = start_std; i < end_std; i++) {

          icelcr[i] = halo->send_list[i] + 1;
          ipercr[i] = t_id;

        } /* End of loop on standard ghost cells */

        if (mesh->halo_type == CS_HALO_EXTENDED) {

          start_ext = halo->perio_lst[shift + 4*rank_id + 2];
          length = halo->perio_lst[shift + 4*rank_id + 3];
          end_ext = start_ext + length;

          for (i = start_ext; i < end_ext; i++) {

            icelcr[i] = halo->send_list[i] + 1;
            ipercr[i] = t_id;

          }

        } /* If an extended halo exists */

      } /* If on local rank */

    } /* End of loop on ranks */

  } /* End of loop on transformations */

}

/*----------------------------------------------------------------------------
 * Apply transformation to the location of a particle.
 *
 * Fortran Interface :
 *
 * SUBROUTINE LAGPER
 * *****************
 *
 * INTEGER          ITRANS        :  -> : transformation id buffer
 * DOUBLE PRECISION VTX_A         :  -> : location of vertex before transform.
 * DOUBLE PRECISION VTX_B         : <-  : location of the vertex after transform.
 *
 * Returns:
 *----------------------------------------------------------------------------*/

void
CS_PROCF (lagper, LAGPER)(const cs_int_t   *itrans,
                          const cs_real_t   vtx_a[],
                                cs_real_t   vtx_b[])
{
  cs_int_t  i, j, rev_id;

  cs_real_t  vect[4];
  cs_real_t  matrix[3][4];

  cs_mesh_t  *mesh = cs_glob_mesh;
  fvm_periodicity_t  *periodicity = mesh->periodicity;

  /* Get the matrix for this transformation */

  rev_id = fvm_periodicity_get_reverse_id(mesh->periodicity, *itrans);
  fvm_periodicity_get_matrix(periodicity, rev_id, matrix);

  /* Initialize vectors */

  for (i = 0; i < 3; i++) {
    vtx_b[i] = 0.;
    vect[i] = vtx_a[i];
  }
  vect[3] = 1;

  /* Compute transformation */

  for (i = 0; i < 3; i++)
    for (j = 0; j < 4; j++)
      vtx_b[i] += matrix[i][j]*vect[j];

}


/*----------------------------------------------------------------------------
 * Apply rotation on the velocity vector of a particle.
 *
 * Fortran Interface :
 *
 * SUBROUTINE LAGVEC
 * *****************
 *
 * INTEGER          ITRANS        :  -> : transformation id
 * DOUBLE PRECISION VECTI         :  -> : vector before transformation
 * DOUBLE PRECISION VECTF         : <-  : vector after transformation
 *
 * Returns:
 *----------------------------------------------------------------------------*/

void
CS_PROCF (lagvec, LAGVEC)(const cs_int_t   *itrans,
                          const cs_real_t  vecti[],
                                cs_real_t  vectf[])
{
  cs_int_t  i, j, rev_id;

  cs_real_t  matrix[3][4];

  cs_mesh_t  *mesh = cs_glob_mesh;
  fvm_periodicity_t  *periodicity = mesh->periodicity;

  /* test if the transformation is a rotation */

  if (FVM_PERIODICITY_ROTATION ==
      fvm_periodicity_get_type(periodicity, *itrans)) {

    /* Get the matrix for this transformation */

    rev_id = fvm_periodicity_get_reverse_id(mesh->periodicity, *itrans);
    fvm_periodicity_get_matrix(periodicity, rev_id, matrix);

    for (i = 0; i < 3; i++) {
      vectf[i] = 0;
      for (j = 0; j < 3; j++)
        vectf[i] += matrix[i][j]*vecti[j];
    }

  } /* If the periodicity is a rotation */

  else {

    /* There is no rotation. Copy the input vector */

    for (i = 0; i < 3; i++)
      vectf[i] = vecti[i];

  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
