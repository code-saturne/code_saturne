#ifndef __CS_ADVECTION_FIELD_H__
#define __CS_ADVECTION_FIELD_H__

/*============================================================================
 * Manage the definition/setting of advection fields
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

#include "cs_cdo_connect.h"
#include "cs_cdo_local.h"
#include "cs_cdo_quantities.h"
#include "cs_mesh_location.h"
#include "cs_param.h"
#include "cs_property.h"
#include "cs_time_step.h"
#include "cs_xdef.h"
#include "cs_xdef_eval.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/


#define CS_ADVECTION_FIELD_POST_COURANT (1 << 0)  // postprocess Courant number
#define CS_ADVECTION_FIELD_STEADY       (1 << 1)  // steady-state field

typedef struct {

  int             id;
  char  *restrict name;

  cs_flag_t   loc_flag;      // Short descriptor for localization
  cs_flag_t   flag;          // Short descriptor dedicated to postprocessing
  int         vtx_field_id;  // id among cs_field_t structures (-1 if not used)
  int         cell_field_id; // id among cs_field_t structures (-1 if not used)

  /* We assume that there is only one definition assicated to an advection
     field */
  cs_xdef_t               *definition;

  /* Function pointers */
  cs_xdef_eval_t          *get_eval_all_vertices;
  cs_xdef_eval_t          *get_eval_at_cell;
  cs_xdef_eval_cw_t       *get_eval_at_cell_cw;
  cs_xdef_eval_cw_xyz_t   *get_eval_at_xyz_cw;

} cs_adv_field_t;

/* List of available keys for setting an advection field */
typedef enum {

  CS_ADVKEY_DEFINE_AT_CELLS,
  CS_ADVKEY_DEFINE_AT_VERTICES,
  CS_ADVKEY_POST_COURANT,
  CS_ADVKEY_STATE_STEADY,
  CS_ADVKEY_N_KEYS

} cs_advection_field_key_t;

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set shared pointers to main domain members
 *
 * \param[in]  quant       additional mesh quantities struct.
 * \param[in]  connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]  time_step   pointer to a time step structure
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_set_shared_pointers(const cs_cdo_quantities_t  *quant,
                                       const cs_cdo_connect_t     *connect,
                                       const cs_time_step_t       *time_step);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the number of allocated cs_adv_field_t structures
 *
 * \return the number of advection fields
 */
/*----------------------------------------------------------------------------*/

int
cs_advection_field_get_n_fields(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Search in the array of advection field structures which one has
 *         the name given in argument
 *
 * \param[in]  name        name of the advection field
 *
 * \return a pointer to a cs_adv_field_t structure or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_adv_field_t *
cs_advection_field_by_name(const char   *name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Search in the array of advection field structures which one has
 *         the id given in argument
 *
 * \param[in]  id        identification number
 *
 * \return a pointer to a cs_adv_field_t structure or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_adv_field_t *
cs_advection_field_by_id(int      id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add and initialize a new advection field structure
 *
 * \param[in]  name        name of the advection field
 *
 * \return a pointer to the new allocated cs_adv_field_t structure
 */
/*----------------------------------------------------------------------------*/

cs_adv_field_t *
cs_advection_field_add(const char   *name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free all alllocated cs_adv_field_t structures and its related array
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_destroy_all(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if the given advection field has the name ref_name
 *
 * \param[in]  adv         pointer to a cs_adv_field_t structure to test
 * \param[in]  ref_name    name of the advection field to find
 *
 * \return true if the name of the advection field is ref_name otherwise false
 */
/*----------------------------------------------------------------------------*/

bool
cs_advection_field_check_name(const cs_adv_field_t   *adv,
                              const char             *ref_name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  returns true if the advection field is uniform, otherwise false
 *
 * \param[in]    adv    pointer to a property to test
 *
 * \return  true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_advection_field_is_uniform(const cs_adv_field_t   *adv);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  returns true if the advection field is uniform in each cell
 *         otherwise false
 *
 * \param[in]    adv    pointer to a property to test
 *
 * \return  true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_advection_field_is_cellwise(const cs_adv_field_t   *adv);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the name of an advection field
 *
 * \param[in]    adv    pointer to an advection field structure
 *
 * \return  the name of the related advection field
 */
/*----------------------------------------------------------------------------*/

const char *
cs_advection_field_get_name(const cs_adv_field_t   *adv);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the type of definition used to set the current advection
 *         field structure
 *
 * \param[in]    adv    pointer to an advection field structure
 *
 * \return  the type of definition
 */
/*----------------------------------------------------------------------------*/

cs_xdef_type_t
cs_advection_field_get_deftype(const cs_adv_field_t   *adv);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Print all setup information related to cs_adv_field_t structures
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_log_setup(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set optional parameters related to a cs_adv_field_t structure
 *
 * \param[in, out]  adv       pointer to a cs_adv_field_t structure
 * \param[in]       key       key related to the member of adv to set
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_set_option(cs_adv_field_t            *adv,
                              cs_advection_field_key_t   key);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the value of a cs_adv_field_t structure
 *
 * \param[in, out]  adv       pointer to a cs_adv_field_t structure
 * \param[in]       vector    values to set
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_def_by_value(cs_adv_field_t    *adv,
                                cs_real_t          vector[3]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_adv_field_t structure thanks to an analytic function
 *
 * \param[in, out]  adv     pointer to a cs_adv_field_t structure
 * \param[in]       func    pointer to a function
 * \param[in]       input   NULL or pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_def_by_analytic(cs_adv_field_t        *adv,
                                   cs_analytic_func_t    *func,
                                   void                  *input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_adv_field_t structure thanks to an array of values
 *
 * \param[in, out]  adv       pointer to a cs_adv_field_t structure
 * \param[in]       loc       information to know where are located values
 * \param[in]       array     pointer to an array
 * \param[in]       index     optional pointer to the array index
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_def_by_array(cs_adv_field_t    *adv,
                                cs_flag_t          loc,
                                cs_real_t         *array,
                                cs_lnum_t         *index);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_adv_field_t structure thanks to an array of values
 *
 * \param[in, out]  adv       pointer to a cs_adv_field_t structure
 * \param[in]       field     pointer to a cs_field_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_def_by_field(cs_adv_field_t    *adv,
                                cs_field_t        *field);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create all needed cs_field_t structures related to an advection
 *         field
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_create_fields(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get a cs_field_t structure related to an advection field and a mesh
 *         location
 *
 * \param[in]  adv         pointer to a cs_adv_field_t structure
 * \param[in]  ml_type     type of mesh location (cells or vertices)
 *
 * \return a pointer to a cs_field_t structure
 */
/*----------------------------------------------------------------------------*/

cs_field_t *
cs_advection_field_get_field(cs_adv_field_t           *adv,
                             cs_mesh_location_type_t   ml_type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of the advection field at the cell center
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      adv     pointer to a cs_adv_field_t structure
 * \param[in, out] vect    pointer to a cs_nvec3_t structure (meas + unitv)
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_in_cell(const cs_cell_mesh_t   *cm,
                           const cs_adv_field_t   *adv,
                           cs_nvec3_t             *vect);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of the advection field at the cell center
 *
 * \param[in]      c_id    id of the current cell
 * \param[in]      adv     pointer to a cs_adv_field_t structure
 * \param[in, out] vect    pointer to a cs_nvec3_t structure (meas + unitv)
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_get_cell_vector(cs_lnum_t               c_id,
                                   const cs_adv_field_t   *adv,
                                   cs_nvec3_t             *vect);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of the advection field for a given face
 *
 * \param[in]      adv     pointer to a cs_adv_field_t structure
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      xyz     coordinates where to evaluate the advection field
 * \param[in, out] vect    pointer to a cs_nvec3_t structure (meas + unitv)
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_get_at_xyz(const cs_adv_field_t   *adv,
                              const cs_cell_mesh_t   *cm,
                              const cs_real_3_t       xyz,
                              cs_nvec3_t             *vect);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of the advection field at cell centers
 *
 * \param[in]      adv           pointer to a cs_adv_field_t structure
 * \param[in, out] cell_values   array of values at cell centers
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_at_cells(const cs_adv_field_t  *adv,
                            cs_real_t             *cell_values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of the advection field at vertices
 *
 * \param[in]      adv          pointer to a cs_adv_field_t structure
 * \param[in, out] vtx_values   array storing the results
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_at_vertices(const cs_adv_field_t  *adv,
                               cs_real_t             *vtx_values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of the flux of the advection field across the
 *         the dual faces of a cell
 *
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in]      adv      pointer to a cs_adv_field_t structure
 * \param[in, out] fluxes   array of values attached to dual faces of a cell
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_get_flux_dfaces(const cs_cell_mesh_t         *cm,
                                   const cs_adv_field_t         *adv,
                                   cs_real_t                    *fluxes);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of the flux of the advection field across the
 *         triangle defined by the two vertices of an edge and the barycenter
 *         of a face.
 *
 * \param[in]  adv       pointer to a cs_adv_field_t structure
 * \param[in]  cm        pointer to a cs_cell_mesh_t structure
 * \param[in]  tef_meas  area of the triangle tef
 * \param[in]  f         id of the face in the current cell
 * \param[in]  e         id of the edge in the current cell
 * \param[in]  v1        id of the first vertex in the current cell
 * \param[in]  v2        id of the second vertex in the current cell
 *
 * \return the value of the flux across tef
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_advection_field_get_flux_tef(const cs_adv_field_t        *adv,
                                const cs_cell_mesh_t        *cm,
                                const cs_real_t              tef_meas,
                                short int                    f,
                                short int                    e,
                                short int                    v1,
                                short int                    v2);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  For each cs_adv_field_t structures, update the values of the related
 *         field(s)
 *
 * \param[in]      cur2prev    true or false
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_update(bool   cur2prev);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the Peclet number in each cell
 *
 * \param[in]      adv        pointer to the advection field struct.
 * \param[in]      diff       pointer to the diffusion property struct.
 * \param[in, out] peclet     pointer to an array storing Peclet number
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_get_peclet(const cs_adv_field_t        *adv,
                        const cs_property_t         *diff,
                        cs_real_t                    peclet[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_ADVECTION_FIELD_H__ */
