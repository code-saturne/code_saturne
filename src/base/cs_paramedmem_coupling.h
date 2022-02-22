#ifndef __CS_PARAMEDMEM_HXX__
#define __CS_PARAMEDMEM_HXX__

/*============================================================================
 * Coupling using ParaMEDMEM
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*============================================================================
 * Structure definitions
 *============================================================================*/

typedef struct _cs_paramedmem_coupling_t cs_paramedmem_coupling_t;

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 *  Global variable definitions
 *============================================================================*/

typedef enum {
  CS_MEDCPL_ON_CELLS,
  CS_MEDCPL_ON_NODES
} cs_medcpl_space_discr_t;

typedef enum {
  CS_MEDCPL_NO_TIME,
  CS_MEDCPL_ONE_TIME,
  CS_MEDCPL_LINEAR_TIME
} cs_medcpl_time_discr_t;

typedef enum {
  CS_MEDCPL_FIELD_INT_CONSERVATION,
  CS_MEDCPL_FIELD_INT_MAXIMUM,
  CS_MEDCPL_FIELD_EXT_CONSERVATION,
  CS_MEDCPL_FIELD_EXT_MAXIMUM,
  CS_MEDCPL_FIELD_N_NATURE
} cs_medcpl_field_nature_t;

/*============================================================================
 * Public C++ function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve coupling struct pointer by id
 *
 * \param[in] cpl_id  index of the sought coupling
 *
 * \return pointer to cs_paramedmem_coupling_t struct. Raise an error if
 * the coupling does not exist.
 */
/*----------------------------------------------------------------------------*/

cs_paramedmem_coupling_t *
cs_paramedmem_coupling_by_id(int cpl_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve coupling struct pointer by name
 *
 * \param[in] name  name of the coupling
 *
 * \return pointer to cs_paramedmem_coupling_t struct or NULL if not found.
 *
 */
/*----------------------------------------------------------------------------*/

cs_paramedmem_coupling_t *
cs_paramedmem_coupling_by_name(const char *name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a new ParaMEDMEM coupling
 *
 * \param[in] app1_name  Name of app n°1 or NULL if calling app is app1
 * \param[in] app2_name  Name of app n°2 or NULL if calling app is app2
 * \param[in] cpl_name   Name of the coupling. If NULL an automatic name is generated.
 *
 * \return pointer to newly created cs_paramedmem_coupling_t structure.
 *
 */
/*----------------------------------------------------------------------------*/

cs_paramedmem_coupling_t *
cs_paramedmem_coupling_create(const char *app1_name,
                              const char *app2_name,
                              const char *cpl_name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy a given ParaMEDMEM coupling structure
 *
 * \param[in] c pointer to cs_paramedmem_coupling_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_paramedmem_coupling_destroy(cs_paramedmem_coupling_t  *c);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy all coupling structures
 */
/*----------------------------------------------------------------------------*/

void
cs_paramedmem_coupling_all_finalize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define coupled mesh based on a selection criteria
 *
 * \param[in] c         pointer to cs_paramedmem_coupling_t struct
 * \param[in] sel_crit  geometrical selection criteria (string)
 * \param[in] elt_dim   dimension of coupled elements
 */
/*----------------------------------------------------------------------------*/

void
cs_paramedmem_add_mesh_from_criteria(cs_paramedmem_coupling_t  *c,
                                     const char                *sel_crit,
                                     int                        elt_dim);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define coupled mesh based on a cs_zone_t pointer
 *
 * \param[in] c     pointer to cs_paramedmem_coupling_t struct
 * \param[in] zone  pointer to cs_zone_t struct
 */
/*----------------------------------------------------------------------------*/

void
cs_paramedmem_add_mesh_from_zone(cs_paramedmem_coupling_t  *c,
                                 const cs_zone_t           *zone);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get number of defined couplings
 *
 * \return number of defined couplings (int)
 */
/*----------------------------------------------------------------------------*/
int
cs_paramedmem_get_number_of_couplings(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get number of elements of coupled mesh
 *
 * \param[in] coupling  pointer to cs_paramedmem_coupling_t struct
 *
 * \return number of elements in mesh associated to coupling
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_paramedmem_mesh_get_n_elts(const cs_paramedmem_coupling_t  *coupling);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get indirection list for elements in coupled mesh
 *
 * \param[in] coupling  pointer to cs_paramedmem_coupling_t struct
 *
 * \return cs_lnum_t pointer to indirection list
 */
/*----------------------------------------------------------------------------*/

const cs_lnum_t *
cs_paramedmem_mesh_get_elt_list(const cs_paramedmem_coupling_t  *coupling);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a coupled field
 *
 * \param[in] c             pointer to cs_paramedmem_coupling_t struct
 * \param[in] name          name of field
 * \param[in] dim           field dimension
 * \param[in] field_nature  field nature flag
 * \param[in] space_discr   field space discretisation (nodes or cells)
 * \param[in] time_discr    field coupling time discretisation
 *
 * \return index of field within the storing vector
 *
 */
/*----------------------------------------------------------------------------*/

int
cs_paramedmem_def_coupled_field(cs_paramedmem_coupling_t  *c,
                                const char                *name,
                                int                        dim,
                                cs_medcpl_field_nature_t   field_nature,
                                cs_medcpl_space_discr_t    space_discr,
                                cs_medcpl_time_discr_t     time_discr);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a coupled field based on a cs_field_t pointer
 *
 * \param[in] c           pointer to cs_paramedmem_coupling_t struct
 * \param[in] f           pointer to cs_field_t struct
 * \param[in] fn          field nature flag
 * \param[in] time_discr  field coupling time discretisation
 *
 * \return index of field within the storing vector
 */
/*----------------------------------------------------------------------------*/

int
cs_paramedmem_def_coupled_field_from_cs_field(cs_paramedmem_coupling_t *c,
                                              cs_field_t               *f,
                                              cs_medcpl_field_nature_t  fn,
                                              cs_medcpl_time_discr_t    td);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Write values before sending operation
 *
 * \param[in] c      pointer to cs_paramedmem_coupling_t structure
 * \param[in] name   name of field
 * \param[in] values array of values to write
 */
/*----------------------------------------------------------------------------*/

void
cs_paramedmem_field_export(cs_paramedmem_coupling_t  *c,
                           const char                *name,
                           const double               values[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read values before sending operation
 *
 * \param[in] c      pointer to cs_paramedmem_coupling_t structure
 * \param[in] name   name of field
 * \param[in] values array in which values will be stored
 */
/*----------------------------------------------------------------------------*/

void
cs_paramedmem_field_import(cs_paramedmem_coupling_t  *c,
                           const char                *name,
                           double                     values[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Sync the coupling's InterpKernelDEC
 *
 * \param[in] c pointer to cs_paramedmem_coupling_t structure
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_paramedmem_sync_dec(cs_paramedmem_coupling_t  *c);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Send values of field attached to DEC
 *
 * \param[in] c pointer to cs_paramedmem_coupling_t structure
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_paramedmem_send_data(cs_paramedmem_coupling_t  *c);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Recieve values of field attached to DEC
 *
 * \param[in] c pointer to cs_paramedmem_coupling_t structure
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_paramedmem_recv_data(cs_paramedmem_coupling_t  *c);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Attach a field to InterpKernelDEC for send operation using its index
 *
 * \param[in] c         pointer to cs_paramedmem_coupling_t structure
 * \param[in] field_id  index of field in storing vector
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_paramedmem_attach_field_by_id(cs_paramedmem_coupling_t  *c,
                                 int                        field_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Attach a field to InterpKernelDEC for send operation using its name
 *
 * \param[in] c     pointer to cs_paramedmem_coupling_t structure
 * \param[in] name  name of field (string)
 */
/*----------------------------------------------------------------------------*/

void
cs_paramedmem_attach_field_by_name(cs_paramedmem_coupling_t  *c,
                                   const char                *name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Send values of a field. If vals pointer is non-null,
 * values are updated before send
 *
 * \param[in] c     pointer to cs_paramedmem_coupling_t structure
 * \param[in] name  name of field
 * \param[in] vals  array of values to write
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_paramedmem_send_field_vals(cs_paramedmem_coupling_t *c,
                              const char               *name,
                              const double             *vals);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Recieve values of a field.
 *
 * \param[in] c     pointer to cs_paramedmem_coupling_t structure
 * \param[in] name  name of field
 * \param[in] vals  array of values to write
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_paramedmem_recv_field_vals(cs_paramedmem_coupling_t *c,
                              const char               *name,
                              double                   *vals);

/*----------------------------------------------------------------------------*/
/*!
 * \brief initialize couplings based on user functions
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_paramedmem_coupling_all_init(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief initialize coupled mesh and fields based on user functions
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_paramedmem_coupling_define_mesh_fields(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log ParaMEDMEM coupling setup information
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_paramedmem_coupling_log_setup(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PARAMEDMEM_HXX__ */
