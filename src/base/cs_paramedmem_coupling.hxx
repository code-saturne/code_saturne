#ifndef __CS_PARAMEDMEM_HXX__
#define __CS_PARAMEDMEM_HXX__

/*============================================================================
 * Coupling using ParaMEDMEM
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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * MED library headers
 *----------------------------------------------------------------------------*/

#if defined(HAVE_PARAMEDMEM)
#include <ParaFIELD.hxx>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

#if defined(HAVE_PARAMEDMEM)
using namespace MEDCoupling;

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*============================================================================
 * Structure definitions
 *============================================================================*/

typedef struct _cs_paramedmem_coupling_t cs_paramedmem_coupling_t;

/*============================================================================
 *  Global variable definitions
 *============================================================================*/

BEGIN_C_DECLS

/*============================================================================
 * Public C++ function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create a paramedmem coupling based on an InterpKernelDEC.
 *
 * The latter is created using the the lists of ranks provided as
 * input to this function.
 *
 * parameters:
 *   name              <-- coupling name
 *   grp1_global_ranks <-- array of ranks of group 1
 *   grp1_size         <-- size of grp1_global_ranks array
 *   grp2_global_ranks <-- array of ranks of group 2
 *   grp2_size         <-- size of grp2_global_ranks array
 *
 * return:
 *   pointer to new coupling object
 *----------------------------------------------------------------------------*/

cs_paramedmem_coupling_t *
cs_paramedmem_interpkernel_create(const char  *name,
                                  int         *grp1_global_ranks,
                                  int          grp1_size,
                                  int         *grp2_global_ranks,
                                  int          grp2_size);

/*----------------------------------------------------------------------------
 * Define new ParaMEDMEM coupling.
 *
 * arguments:
 *   name     <-- name of coupling
 *   send_dec <-- send Data Exchange Channel
 *   recv_dec <-- receive Data Exchange Channel
 *----------------------------------------------------------------------------*/

void
cs_paramedmem_destroy(cs_paramedmem_coupling_t  **coupling);

/*----------------------------------------------------------------------------
 * Define nodal mesh for ParaMEDMEM coupling from selection criteria.
 *
 * parameters:
 *   coupling        <-- partially initialized ParaMEDMEM coupling structure
 *   name            <-- name of coupling mesh
 *   select_criteria <-- selection criteria
 *   elt_dim         <-- element dimension
 *   is_source       <-- true if fields located on mesh are sent
 *   is_dest         <-- true if fields located on mesh are received
 *
 * returns:
 *   id of created mesh in coupling
 *----------------------------------------------------------------------------*/

int
cs_paramedmem_define_mesh(cs_paramedmem_coupling_t  *coupling,
                          const char                *name,
                          const char                *select_criteria,
                          int                        elt_dim,
                          bool                       is_source,
                          bool                       is_dest);

/*----------------------------------------------------------------------------
 * Initialize nodal coupled meshes.
 *
 * parameters:
 *   coupling        <-- partially initialized ParaMEDMEM coupling structure
 *----------------------------------------------------------------------------*/

void
cs_paramedmem_init_meshes(cs_paramedmem_coupling_t  *coupling);

/*----------------------------------------------------------------------------
 * Return the ParaMEDMEM mesh id associated with a given mesh name,
 * or -1 if no association found.
 *
 * parameters:
 *   coupling  <-- coupling structure
 *   mesh_name <-- mesh name
 *
 * returns:
 *    mesh id for this coupling, or -1 if mesh name is not associated
 *    with this coupling.
 *----------------------------------------------------------------------------*/

int
cs_paramedmem_mesh_id(cs_paramedmem_coupling_t  *coupling,
                      const char                *mesh_name);

/*----------------------------------------------------------------------------
 * Get number of associated coupled elements in coupled mesh
 *
 * parameters:
 *   coupling <-- ParaMEDMEM coupling structure
 *   mesh_id  <-- id of coupled mesh in coupling
 *
 * returns:
 *   number of elements in coupled mesh
 *----------------------------------------------------------------------------*/

cs_lnum_t
cs_paramedmem_mesh_get_n_elts(const cs_paramedmem_coupling_t *coupling,
                              int                             mesh_id);

/*----------------------------------------------------------------------------
 * Get local list of coupled elements (0 to n-1 numbering) for a coupled mesh
 *
 * parameters:
 *   coupling <-- ParaMEDMEM coupling structure
 *   mesh_id  <-- id of coupled mesh in coupling
 *----------------------------------------------------------------------------*/

const cs_lnum_t *
cs_paramedmem_mesh_get_elt_list(const cs_paramedmem_coupling_t *coupling,
                                int                             mesh_id);

/*----------------------------------------------------------------------------
 * Create a MEDCoupling field structure.
 *
 * parameters:
 *   coupling  <-- MED coupling structure.
 *   name      <-- field name.
 *   mesh_id   <-- id of associated mesh in structure.
 *   dim       <-- number of field components.
 *   type      <-- mesh mesh (ON_NODES, ON_CELLS)
 *   td        <-- time discretization type
 *   dirflag   <-- 1: send, 2: receive
 *
 * returns
 *   field id in coupling structure
 *----------------------------------------------------------------------------*/

int
cs_paramedmem_field_add(cs_paramedmem_coupling_t  *coupling,
                        const char                *name,
                        int                        mesh_id,
                        int                        dim,
                        TypeOfField                type,
                        TypeOfTimeDiscretization   td,
                        int                        dirflag);

/*----------------------------------------------------------------------------
 * Return the ParaMEDMEM field id associated with given mesh and field names,
 * or -1 if no association found.
 *
 * parameters:
 *   coupling <-- coupling structure.
 *   mesh_id  <-- id of associated mesh in structure.
 *   name     <-- field name.
 *
 * returns
 *   field id in coupling structure, or -1 if not found
 *----------------------------------------------------------------------------*/

int
cs_paramedmem_field_get_id(cs_paramedmem_coupling_t  *coupling,
                           int                        mesh_id,
                           const char                *name);

/*----------------------------------------------------------------------------
 * Return ParaMEDMEM::ParaFIELD object associated with a given field id.
 *
 * parameters:
 *   coupling  <-- pointer to associated coupling
 *   field_id  <-- id of associated field structure
 *
 * returns:
 *   pointer to ParaFIELD to which values were assigned
 *----------------------------------------------------------------------------*/

MEDCoupling::ParaFIELD *
cs_paramedmem_field_get(cs_paramedmem_coupling_t  *coupling,
                        int                        field_id);

/*----------------------------------------------------------------------------
 * Write field associated with a mesh to MEDCoupling.
 *
 * Assigning a negative value to the time step indicates a time-independent
 * field (in which case the time_value argument is unused).
 *
 * parameters:
 *   coupling     <-- pointer to associated coupling
 *   field_id     <-- id of associated field
 *   on_parent    <-- if true, values are defined on parent mesh
 *   field_values <-- array of associated field value arrays
 *----------------------------------------------------------------------------*/

void
cs_paramedmem_field_export(cs_paramedmem_coupling_t  *coupling,
                           int                        field_id,
                           bool                       on_parent,
                           const double               field_values[]);

/*----------------------------------------------------------------------------
 * Read field associated with a mesh from MEDCoupling.
 *
 * Only double precision floating point values are considered.
 *
 * Assigning a negative value to the time step indicates a time-independent
 * field (in which case the time_value argument is unused).
 *
 * parameters:
 *   coupling     <-- pointer to associated coupling
 *   field_id     <-- id of associated field
 *   on_parent    <-- if true, values are defined on parent mesh
 *   field_values <-- array of associated field value arrays
 *----------------------------------------------------------------------------*/

void
cs_paramedmem_field_import(cs_paramedmem_coupling_t  *coupling,
                           int                        field_id,
                           bool                       on_parent,
                           double                     field_values[]);

/*----------------------------------------------------------------------------
 * Synchronize DEC assciated with a given coupling.
 *
 * This sync function needs to be called at least once before exchanging data.
 * dec->synchronize() creates the interpolation matrix between the two codes!
 *
 * parameters:
 *   coupling    <-- coupling structure.
 *   dec_to_sync <-- 1 for send_dec, != 1 for recv_dec
 *----------------------------------------------------------------------------*/

void
cs_paramedmem_sync_dec(cs_paramedmem_coupling_t  *coupling,
                       int                        dec_to_sync);

/*----------------------------------------------------------------------------
 * Send the values related to a coupling
 *
 * parameters:
 *   coupling <-> coupling structure.
 *----------------------------------------------------------------------------*/

void
cs_paramedmem_send_data(cs_paramedmem_coupling_t  *coupling);

/*----------------------------------------------------------------------------
 * Receive the values related to a coupling
 *
 * parameters:
 *   coupling <-> coupling structure.
 *----------------------------------------------------------------------------*/

void
cs_paramedmem_recv_data(cs_paramedmem_coupling_t  *coupling);

/*----------------------------------------------------------------------------
 * Link a given field to the DEC before send/recv
 *
 * parameters:
 *   coupling <-> coupling structure.
 *   field_id <-> associated field id
 *----------------------------------------------------------------------------*/

void
cs_paramedmem_reattach_field(cs_paramedmem_coupling_t  *coupling,
                             int                        field_id);

/*============================================================================
 * Public C++ function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Map MPI ranks within cs_glob_mpi_comm to their values in MPI_COMM_WORLD.
 *
 * The caller is responsible for freeing the returned array
 *
 * return:
 *   list of ranks in MPI_COMM_WORLD
 *----------------------------------------------------------------------------*/

int *
cs_paramedmem_get_mpi_comm_world_ranks(void);

/*----------------------------------------------------------------------------*/

#endif

END_C_DECLS

#endif /* __CS_PARAMEDMEM_HXX__ */
