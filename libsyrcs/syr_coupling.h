/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2010 EDF S.A., France
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

#ifndef _SYR_COUPLING_H_
#define _SYR_COUPLING_H_

/*============================================================================
 * Main API functions for coupling between Syrthes and Code_Saturne
 *
 * Library: Code_Saturne                               Copyright EDF 2006-2008
 *============================================================================*/

/*----------------------------------------------------------------------------
 * System headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "syr_comm.h"

/*----------------------------------------------------------------------------
 * Structure definitions
 *----------------------------------------------------------------------------*/

typedef struct _syr_coupling_t syr_coupling_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize syr_coupling_t structure
 *
 * arguments:
 *   coupling_id <-- Id of Syrthes coupling (0 to n-1)
 *   cs_app_name <-- Application name of Code_Saturne MPI process
 *   sock_str    <-- hostname:socknum of first coupled
 *                   Code_Saturne process, or NULL
 *   comm_type   <-- Type of comunication used
 *   comm_echo   <-- Optional echo to standard output
 *----------------------------------------------------------------------------*/

syr_coupling_t  *
syr_coupling_initialize(int               coupling_id,
                        const char       *cs_app_name,
                        const char       *sock_str,
                        syr_comm_type_t   comm_type,
                        int               comm_echo);

/*----------------------------------------------------------------------------
 * Finalize syr_coupling_t structure
 *----------------------------------------------------------------------------*/

syr_coupling_t  *
syr_coupling_finalize(syr_coupling_t  *coupling);

/*----------------------------------------------------------------------------
 * Receive Code_Saturne's coupled surface mesh
 *
 * the elt_connect and vtx_coords arrays are allocated here.
 *
 * arguments:
 *   coupling    <-- Associated coupling object
 *   sp_dim      --> Spacial dimension for fluid
 *   n_vtx,      --> Number of vertices
 *   n_elts,     --> Number of segments
 *   vtx_coords  --> Vertex coordinates
 *   elt_connect --> Segment or triangle connectivity
 *----------------------------------------------------------------------------*/

void
syr_coupling_receive_bc_mesh(syr_coupling_t  *coupling,
                             int             *sp_dim,
                             int             *n_vtx,
                             int             *n_elts,
                             double         **vtx_coords,
                             int            **elt_connect);

/*----------------------------------------------------------------------------
 * Exchange of synchronization (supervision) messages
 *
 * parameters:
 *  coupling <-- Associated coupling object
 *  is_last  --> Last time step or iteration indicator
 *  is_end   --> Calculation stop indicator
 *----------------------------------------------------------------------------*/

void
syr_coupling_supervise(syr_coupling_t  *coupling,
                       int             *is_last,
                       int             *is_end);

/*----------------------------------------------------------------------------
 * Data exchange prior to iteration
 *
 * Send wall temperature
 * Receive fluid temperature and pseudo-exchange coefficient
 * Possibly receive CFD code time step
 *
 * parameters:
 *   coupling <-- Associated coupling object
 *   tpf      <-> Wall Temperature in, Fluid temperature out
 *   hht      --> Pseudo-exchange coefficient
 *   dtfluid  --> Time step set by CFD code if received, -1 otherwise
 *----------------------------------------------------------------------------*/

void
syr_coupling_exchange_var(syr_coupling_t  *coupling,
                          double          *tpf,
                          double          *hht,
                          double          *dtfluid);

/*----------------------------------------------------------------------------*/

#endif /* _SYR_COUPLING_H_ */
