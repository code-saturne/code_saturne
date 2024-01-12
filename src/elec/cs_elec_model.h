#ifndef __CS_ELEC_MODEL_H__
#define __CS_ELEC_MODEL_H__

/*============================================================================
 * General parameters management.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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

#include "cs_defs.h"

#include "cs_domain.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Structure to read properties in dp_ELE
 *----------------------------------------------------------------------------*/

typedef struct {

  int         ngaz;
  int         npoint;
  cs_real_t  *th;
  cs_real_t  *ehgaz;
  cs_real_t  *rhoel;
  cs_real_t  *cpel;
  cs_real_t  *sigel;
  cs_real_t  *visel;
  cs_real_t  *xlabel;
  cs_real_t  *xkabel;

} cs_data_elec_t;

/*----------------------------------------------------------------------------
 * Structure to read transformer parameters in dp_ELE
 *----------------------------------------------------------------------------*/

typedef struct {

  int         nbelec;
  int        *ielecc;
  int        *ielect;
  int        *ielecb;
  int         nbtrf;
  int         ntfref;
  int        *ibrpr;
  int        *ibrsec;
  cs_real_t  *tenspr;
  cs_real_t  *rnbs;
  cs_real_t  *zr;
  cs_real_t  *zi;
  cs_real_t  *uroff;
  cs_real_t  *uioff;

} cs_data_joule_effect_t;

/*----------------------------------------------------------------------------
 * Electrical model options descriptor
 *----------------------------------------------------------------------------*/

typedef struct {

  int         ixkabe;
  int         ntdcla;
  int         irestrike;
  cs_real_t   restrike_point[3];
  cs_real_t   crit_reca[5];
  int         ielcor;
  int         modrec;
  int         idreca;
  int        *izreca;
  cs_real_t   couimp;
  cs_real_t   pot_diff;
  cs_real_t   puisim;
  cs_real_t   coejou;
  cs_real_t   elcou;
  cs_real_t   srrom;

} cs_elec_option_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Pointer to electrical model options structure */

extern const cs_elec_option_t        *cs_glob_elec_option;
extern const cs_data_elec_t          *cs_glob_elec_properties;
extern const cs_data_joule_effect_t  *cs_glob_transformer;

/* Constant for electrical models */

extern const cs_real_t cs_elec_permvi;
extern const cs_real_t cs_elec_epszer;

/*=============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

void
CS_PROCF (elini1, ELINI1) (void);

void
CS_PROCF (elflux, ELFLUX) (int      *iappel);

void
CS_PROCF (elthht, ELTHHT) (int       *mode,
                           cs_real_t *ym,
                           cs_real_t *enthal,
                           cs_real_t *temp);

void
CS_PROCF (ellecd, ELLECD) (void);

void
CS_PROCF (elphyv, ELPHYV) (void);

void
CS_PROCF (eltssc, ELTSSC) (const int       *isca,
                           cs_real_t       *smbrs);

void
CS_PROCF (eliniv, ELINIV) (int      *isuite);

void
CS_PROCF (elreca, ELRECA) (cs_real_t *dt);

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Provide access to cs_elec_option
 *----------------------------------------------------------------------------*/

cs_elec_option_t *
cs_get_glob_elec_option(void);

/*----------------------------------------------------------------------------
 * Provide access to cs_glob_transformer
 *----------------------------------------------------------------------------*/

cs_data_joule_effect_t *
cs_get_glob_transformer(void);

/*----------------------------------------------------------------------------
 * Initialize structures for electrical model
 *----------------------------------------------------------------------------*/

void
cs_electrical_model_initialize(void);

/*----------------------------------------------------------------------------
 * Destroy structures for electrical model
 *----------------------------------------------------------------------------*/

void
cs_electrical_model_finalize(void);

/*----------------------------------------------------------------------------
 * Specific initialization for electric arc
 *----------------------------------------------------------------------------*/

void
cs_electrical_model_specific_initialization(void);

/*----------------------------------------------------------------------------
 * Read properties file
 *----------------------------------------------------------------------------*/

void
cs_electrical_properties_read(void);

/*----------------------------------------------------------------------------
 * compute specific electric arc fields
 *----------------------------------------------------------------------------*/

void
cs_elec_compute_fields(const cs_mesh_t  *mesh,
                       int               call_id);

/*----------------------------------------------------------------------------
 * compute physical properties
 *----------------------------------------------------------------------------*/

void
cs_elec_physical_properties(cs_domain_t             *domain);

/*----------------------------------------------------------------------------
 * compute source terms for energy
 *----------------------------------------------------------------------------*/

void
cs_elec_source_terms(const cs_mesh_t             *mesh,
                     const cs_mesh_quantities_t  *mesh_quantities,
                     int                          f_id,
                     cs_real_t                   *smbrs);

/*----------------------------------------------------------------------------
 * compute source terms for vector potential
 *----------------------------------------------------------------------------*/

void
cs_elec_source_terms_v(const cs_mesh_t             *mesh,
                       const cs_mesh_quantities_t  *mesh_quantities,
                       int                          f_id,
                       cs_real_3_t                 *smbrv);

/*----------------------------------------------------------------------------
 * add variables fields
 *----------------------------------------------------------------------------*/

void
cs_elec_add_variable_fields(void);

/*----------------------------------------------------------------------------
 * add properties fields
 *----------------------------------------------------------------------------*/

void
cs_elec_add_property_fields(void);

/*----------------------------------------------------------------------------
 * initialize electric fields
 *----------------------------------------------------------------------------*/

void
cs_elec_fields_initialize(const cs_mesh_t  *mesh,
                          int               isuite);

/*----------------------------------------------------------------------------
 * scaling electric quantities
 *----------------------------------------------------------------------------*/

void
cs_elec_scaling_function(const cs_mesh_t             *mesh,
                         const cs_mesh_quantities_t  *mesh_quantities,
                         cs_real_t                   *dt);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Convert enthalpy to temperature at all boundary faces.
 *
 * This handles both user and model enthalpy conversions, so can be used
 * safely whenever conversion is needed.
 *
 * \param[in]   h   enthalpy values
 * \param[out]  t   temperature values
 */
/*----------------------------------------------------------------------------*/

void
cs_elec_convert_h_to_t_faces(const cs_real_t  h[],
                             cs_real_t        t[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Convert single enthalpy value to temperature.
 *
 * \param[in]       ym      mass fraction for each gas
 * \param[in, out]  enthal  enthlapy value
 *
 * \return  temperature value
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_elec_convert_h_to_t(const cs_real_t  ym[],
                       cs_real_t        enthal);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Convert temperature to enthalpy at all cells
 *
 * This handles both user and model temperature conversions, so can be used
 * safely whenever conversion is needed.
 *
 * \param[in]   t   temperature values
 * \param[out]  h   enthalpy values
 */
/*----------------------------------------------------------------------------*/

void
cs_elec_convert_t_to_h_cells(const cs_real_t  t[],
                             cs_real_t        h[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Convert temperature to enthalpy at selected boundary faces.
 *
 * This handles both user and model temperature conversions, so can be used
 * safely whenever conversion is needed.
 *
 * \param[in]   n_faces   number of selected faces
 * \param[in]   face_ids  ids of selected faces
 * \param[in]   t         temperature values (defined on all boundary faces)
 * \param[out]  h         enthalpy values (defined on all boundary faces)
 */
/*----------------------------------------------------------------------------*/

void
cs_elec_convert_t_to_h_faces(const cs_lnum_t  n_faces,
                             const cs_lnum_t  face_ids[],
                             const cs_real_t  t[],
                             cs_real_t        h[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Convert single temperature value to enthalpy.
 *
 * \param[in]       ym    mass fraction for each gas
 * \param[in, out]  temp  temperature value
 *
 * \return  enthalpy values
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_elec_convert_t_to_h(const cs_real_t ym[],
                       cs_real_t       temp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create or access function objects specific to
 *        electric arcs models.
 */
/*----------------------------------------------------------------------------*/

void
cs_elec_define_functions(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_ELEC_MODEL_H__ */
