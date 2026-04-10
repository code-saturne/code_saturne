#ifndef __CS_PARAMEDMEM_HXX__
#define __CS_PARAMEDMEM_HXX__

/*============================================================================
 * MEDCoupling ParaMESH/ParaFIELD wrapper functions.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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

#include "base/cs_defs.h"

#include "base/cs_field.h"
#include "base/cs_zone.h"

/*============================================================================
 * Structure definitions
 *============================================================================*/

struct _cs_paramedmem_coupling_t;
typedef struct _cs_paramedmem_coupling_t cs_paramedmem_coupling_t;

/*============================================================================
 *  Global variable definitions
 *============================================================================*/

typedef enum {
  CS_MEDCPL_ON_CELLS,
  CS_MEDCPL_ON_NODES,
  CS_MEDCPL_ON_NODES_FE
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

typedef enum {
  CS_MEDCPL_INTERPKERNELDEC,
  CS_MEDCPL_CFEMDEC,
} cs_medcpl_dec_t;

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
 * \param[in] cpl_name   Name of the coupling.
 *                       If NULL an automatic name is generated.
 * \param[in] type_dec   Type of DEC used for Data exchange
 *
 * \return pointer to newly created cs_paramedmem_coupling_t structure.
 */
/*----------------------------------------------------------------------------*/

cs_paramedmem_coupling_t *
cs_paramedmem_coupling_create(const char     *app1_name,
                              const char     *app2_name,
                              const char     *cpl_name,
                              cs_medcpl_dec_t type_dec);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a new ParaMEDMEM handler structure with no actual coupling.
 *
 * This can be useful for a "dry run" when setting up a coupling, so as to
 * first debug local commands before actually running in coupled mode.
 *
 * In this case, data "received" matches the initialized values.
 *
 * \param[in] cpl_name   Name of the coupling.
 *                       If NULL an automatic name is generated.
 *
 * \return pointer to newly created cs_paramedmem_coupling_t structure.
 */
/*----------------------------------------------------------------------------*/

cs_paramedmem_coupling_t *
cs_paramedmem_coupling_create_uncoupled(const char *cpl_name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy a given ParaMEDMEM coupling structure
 *
 * \param[in] c pointer to cs_paramedmem_coupling_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_paramedmem_coupling_destroy(cs_paramedmem_coupling_t *c);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy all coupling structures
 */
/*----------------------------------------------------------------------------*/

void
cs_paramedmem_coupling_all_finalize(void);

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
 * \brief initialize couplings based on user functions
 */
/*----------------------------------------------------------------------------*/

void
cs_paramedmem_coupling_all_init(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief initialize coupled mesh and fields based on user functions
 */
/*----------------------------------------------------------------------------*/

void
cs_paramedmem_coupling_define_mesh_fields(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log ParaMEDMEM coupling setup information
 */
/*----------------------------------------------------------------------------*/

void
cs_paramedmem_coupling_log_setup(void);

#endif /* __CS_PARAMEDMEM_HXX__ */
