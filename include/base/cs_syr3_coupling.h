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

#ifndef __CS_SYR_COUPLING_H__
#define __CS_SYR_COUPLING_H__

/*============================================================================
 * Syrthes coupling
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_syr3_comm.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*============================================================================
 * Structure definition
 *============================================================================*/

/* Structure associated to Syrthes coupling */

typedef struct _cs_syr_coupling_t  cs_syr_coupling_t;

/*============================================================================
 *  Global variables definition
 *============================================================================*/

/*============================================================================
 *  Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Get number of Syrthes couplings.
 *
 * Fortran Interface:
 *
 * SUBROUTINE NBCSYR
 * *****************
 *
 * INTEGER          n_couplings     : <-- : number of Syrthes couplings
 *----------------------------------------------------------------------------*/

void CS_PROCF(nbcsyr, NBCSYR)
(
 cs_int_t  *const n_couplings
);

/*----------------------------------------------------------------------------
 * Create nodal coupled mesh.
 * Send vertices's coordinates and connectivity of coupled mesh.
 *
 * Fortran Interface:
 *
 * SUBROUTINE GEOSYR
 * *****************
 *
 * INTEGER          n_couplings     : <-- : number of Syrthes couplings
 *----------------------------------------------------------------------------*/

void CS_PROCF(geosyr, GEOSYR)
(
 cs_int_t  *const n_couplings
);

/*----------------------------------------------------------------------------
 * Get number of boundary faces coupled with Syrthes.
 *
 * Fortran Interface:
 *
 * SUBROUTINE NBFSYR
 * *****************
 *
 * INTEGER          coupl_num       : --> : coupling number
 * INTEGER          n_coupl_faces   : <-- : number of coupled boundary faces
 *----------------------------------------------------------------------------*/

void CS_PROCF(nbfsyr, NBFSYR)
(
 const cs_int_t  *const coupl_num,
       cs_int_t  *const n_coupl_faces
);

/*----------------------------------------------------------------------------
 * Get local numbering of coupled faces
 *
 * Fortran interface:
 *
 * SUBROUTINE LFASYR
 * *****************
 *
 * INTEGER          coupl_num       : --> : coupling number
 * INTEGER          n_coupl_faces   : --> : number of coupled boundary faces
 * INTEGER          coupl_face_list : <-- : list of coupled boundary faces
 *----------------------------------------------------------------------------*/

void CS_PROCF(lfasyr, LFASYR)
(
 const cs_int_t  *const coupl_num,
 const cs_int_t  *const n_coupl_faces,
       cs_int_t  *const coupl_face_list
);

/*----------------------------------------------------------------------------
 * Initialize post processing of Syrthes couplings.
 *
 * Fortran Interface:
 *
 * SUBROUTINE PSTISY
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF(pstisy, PSTISY)
(
 void
);

/*----------------------------------------------------------------------------
 * Get the local (negative) numbers associated with the first and last
 * post processing meshes dedicated to Syrthes couplings
 *
 * Fortran interface:
 *
 * SUBROUTINE PSTESY
 * *****************
 *
 * INTEGER          first_id        : <-- : id of first post processing mesh
 * INTEGER          last_id         : <-- : id of last post processing mesh
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstesy, PSTESY)
(
 cs_int_t  *const first_id,
 cs_int_t  *const last_id
);

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Get number of Syrthes couplings.
 *
 * returns:
 *   number of Syrthes couplings
 *----------------------------------------------------------------------------*/

cs_int_t
cs_syr_coupling_n_couplings(void);

/*----------------------------------------------------------------------------
 * Get pointer to Syrthes coupling.
 *
 * parameters:
 *   coupling_id  -->  Id (0 to n-1) of Syrthes coupling
 *
 * returns:
 *   pointer to Syrthes coupling structure
 *----------------------------------------------------------------------------*/

cs_syr_coupling_t *
cs_syr_coupling_by_id(cs_int_t coupling_id);

/*----------------------------------------------------------------------------
 * Get communicator type associated with Syrthes coupling
 *
 * parameters:
 *   syr_coupling        -->  Syrthes coupling structure
 *
 * returns:
 *   communicator type
 *----------------------------------------------------------------------------*/

cs_comm_type_t
cs_syr_coupling_get_comm_type(const cs_syr_coupling_t *syr_coupling);

/*----------------------------------------------------------------------------
 * Get sending communicator associated with Syrthes coupling
 *
 * parameters:
 *   syr_coupling        -->  coupling structure with Syrthes
 *
 * returns:
 *   pointer to send communicator
 *----------------------------------------------------------------------------*/

cs_comm_t *
cs_syr_coupling_get_send_comm(const cs_syr_coupling_t *syr_coupling);

/*----------------------------------------------------------------------------
 * Get receiving communicator associated with Syrthes coupling
 *
 * parameters:
 *   syr_coupling        -->  coupling structure with Syrthes
 *
 * returns:
 *   pointer to receive communicator
 *----------------------------------------------------------------------------*/

cs_comm_t *
cs_syr_coupling_get_recv_comm(const cs_syr_coupling_t *syr_coupling);

/*----------------------------------------------------------------------------
 * Get number of vertices in coupled mesh
 *
 * parameters:
 *   syr_coupling        -->  Syrthes coupling structure
 *
 * returns:
 *   number of vertices in coupled mesh
 *----------------------------------------------------------------------------*/

cs_int_t
cs_syr_coupling_get_n_vertices(const cs_syr_coupling_t *syr_coupling);

/*----------------------------------------------------------------------------
 * Create a syr_coupling_t structure
 *
 * parameters:
 *   dim                 -->  spatial mesh dimension
 *   ref_axis            -->  reference axis
 *   inv_sel             -->  invert selected faces or not
 *   n_colors            -->  number of colors
 *   colors              -->  color list
 *   n_groups            -->  number of groups
 *   groups              -->  group list
 *   syr_proc_rank       -->  syrthes processus rank
 *   comm_type           -->  communicator type
 *----------------------------------------------------------------------------*/

void
cs_syr_coupling_add(cs_int_t       dim,
                    cs_int_t       ref_axis,
                    cs_bool_t      invsel,
                    cs_int_t       n_colors,
                    cs_int_t      *colors,
                    cs_int_t       n_groups,
                    char         **groups,
#if defined (_CS_HAVE_MPI)
                    cs_int_t       syr_proc_rank,
#endif
                    cs_comm_type_t comm_type);

/*----------------------------------------------------------------------------
 * Initialize communicator for Syrthes coupling
 *
 * parameters:
 *   syr_coupling        -->  Syrthes coupling structure
 *   num_syr_coupling    -->  syrthes coupling number
 *   comm_echo           -->  optional echo to standard output
 *----------------------------------------------------------------------------*/

void
cs_syr_coupling_init_comm(cs_syr_coupling_t *syr_coupling,
                          cs_int_t           num_syr_coupling,
                          cs_int_t           comm_echo);

/*----------------------------------------------------------------------------
 * Destroy cs_syr_coupling_t structures
 *----------------------------------------------------------------------------*/

void
cs_syr_coupling_all_destroy(void);

/*----------------------------------------------------------------------------
 * Define coupled mesh and send it to Syrthes
 *
 * parameters:
 *   syr_coupling        -->  Syrthes coupling structure
 *   coupl_num           -->  syrthes coupling number
 *----------------------------------------------------------------------------*/

void
cs_syr_coupling_init_mesh(cs_syr_coupling_t *syr_coupling,
                          const cs_int_t     coupl_num);

/*----------------------------------------------------------------------------
 * Interpolate a vertex field to an element-centered field
 *
 * parameters:
 *   syr_coupling        -->  Syrthes coupling structure
 *   vtx_values          -->  values defined on vertices
 *   elt_values          <->  values defined on elements
 *----------------------------------------------------------------------------*/

void
cs_syr_coupling_vtx_to_elt(cs_syr_coupling_t        *syr_coupling,
                           cs_real_t          *const vtx_values,
                           cs_real_t                *elt_values);

/*----------------------------------------------------------------------------
 * Interpolate an element-centered field to a vertex field.
 *
 * The size of vtx_values array must be twice the number of vertices.
 * The first half gets values and the second half is used as a working array.
 * The two parts must be contiguous in parallel mode for MPI transfers.
 *
 * parameters:
 *   syr_coupling        -->  Syrthes coupling structure
 *   elt_values          <->  array of values defined on elements
 *   n_vtx_values        -->  number of values defined on vertices
 *   vtx_values          <->  array of values defined on vertices
 *----------------------------------------------------------------------------*/

void
cs_syr_coupling_elt_to_vtx(cs_syr_coupling_t        *syr_coupling,
                           cs_real_t          *const elt_values,
                           cs_int_t                  n_vertices,
                           cs_real_t                *vtx_values);

/*----------------------------------------------------------------------------
 * Update post-processing variables of a Syrthes coupling
 *
 * parameters:
 *   syr_coupling        -->  Syrthes coupling structure
 *   step                -->  0: var = wall temperature
 *                            1: var = fluid temperature
 *                            2: var = exchange coefficient
 *   var                 -->  Pointer to variable values
 *----------------------------------------------------------------------------*/

void
cs_syr_coupling_post_var_update(cs_syr_coupling_t *syr_coupling,
                                int                step,
                                const cs_real_t   *var);

/*----------------------------------------------------------------------------
 * Get the local (negative) numbers associated with the first and last
 * post processing meshes dedicated to Syrthes couplings
 *
 * parameters:
 *   first_mesh_id       <--  Id of first post processing mesh
 *   last_mesh_id        <--  Id of last post processing mesh
 *----------------------------------------------------------------------------*/

void
cs_syr_coupling_post_id_extents(cs_int_t  *const id_mesh_start,
                                cs_int_t  *const id_mesh_end);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CS_SYR_COUPLING_H__ */
