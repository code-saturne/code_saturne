#ifndef __CS_MEG_PROTOTYPES_H__
#define __CS_MEG_PROTOTYPES_H__

/*============================================================================
 * Prototypes for MEG (Mathematical Expression Generator) functions
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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_base.h"
#include "base/cs_field.h"
#include "base/cs_volume_zone.h"
#include "base/cs_boundary_zone.h"
#include "base/cs_ibm.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

typedef void
(cs_ibm_volume_func_t)(const cs_lnum_t    c_id,
                       const cs_real_3_t  xyz,
                       const cs_real_t    t,
                       cs_real_t         *retval);

typedef void
(cs_ibm_fsi_func_t)(cs_real_t  *retval);

/*============================================================================
 * MEG function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function is used to compute user defined values for fields over
 *        a given boundary zone. The mathematical expression is defined in the
 *        GUI.
 */
/*----------------------------------------------------------------------------*/

void
cs_meg_boundary_function
(
  const char       *zone_name,  /*!<[in]  name of a boundary zone */
  const cs_lnum_t   n_elts,     /*!<[in]  number of elements */
  const cs_lnum_t  *elt_ids,    /*!<[in]  list of element ids */
  const cs_real_t   xyz[][3],   /*!<[in]  list of coordinates */
  const char       *field_name, /*!<[in]  name of the variable field */
  const char       *condition,  /*!<[in]  condition type defined as a string */
  cs_real_t        *retvals     /*!<[out] array of computed values */
);

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function is used to compute user defined values for fields over
 *        a given volume zone. The mathematical expression is defined in the
 *        GUI.
 */
/*----------------------------------------------------------------------------*/

void
cs_meg_volume_function
(
  const char      *zone_name,    /*!<[in] name of a volume zone */
  const cs_lnum_t  n_elts,       /*!<[in] number of elements */
  const cs_lnum_t *elt_ids,      /*!<[in] list of element ids */
  const cs_real_t  xyz[][3],     /*!<[in] list of coordinates */
  const char      *fields_names, /*!<[in]  associated variable field names */
  cs_real_t       *fvals[]       /*!<[in, out] array of value arrays */
);

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function is used for the initalization of fields over a given
 *        volume zone. The mathematical expression is defined in the GUI.
 *
 * The caller is responsible for freeing the associated array.
 */
/*----------------------------------------------------------------------------*/

void
cs_meg_initialization
(
  const char      *zone_name,  /*!<[in]  name of a volume zone */
  const cs_lnum_t  n_elts,     /*!<[in]  number of elements */
  const cs_lnum_t *elt_ids,    /*!<[in]  list of element ids */
  const cs_real_t  xyz[][3],   /*!<[in]  list of coordinates */
  const char      *field_name, /*!<[in]  associated variable field name */
  cs_real_t       *retvals     /*!<[out] array of computed values */
);

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function is used to compute source terms over a volume zone. The
 *        mathematical expression is defined in the GUI.
 *
 * The caller is responsible for freeing the returned array.
 */
/*----------------------------------------------------------------------------*/

void
cs_meg_source_terms
(
  const char      *zone_name,   /*!<[in]  name of a volume zone */
  const cs_lnum_t  n_elts,      /*!<[in]  number of elements */
  const cs_lnum_t *elt_ids,     /*!<[in]  list of element ids */
  const cs_real_t  xyz[][3],    /*!<[in]  list of coordinates */
  const char      *name,        /*!<[in]  variable field name */
  const char      *source_type, /*!<[in]  source term type */
  cs_real_t       *retvals      /*!<[out] array of computed values */
);

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function is used to query FSI internal coupling structure values
 *        for a given boundary and structure.
 */
/*----------------------------------------------------------------------------*/

void
cs_meg_fsi_struct
(
  const char       *object_type, /*!<[in] name of object type */
  const char       *name,        /*!<[in] name of matching boundary */
  const cs_real_t   fluid_f[],   /*!<[in] array of fluid forces on the object */
  cs_real_t         val[]        /*!<[in,out] matrix or vector coefficients */
);

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function is used to activate postprocessing writers.
 */
/*----------------------------------------------------------------------------*/

void
cs_meg_post_activate(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function is used to define profile coordinates.
 */
/*----------------------------------------------------------------------------*/

void
cs_meg_post_profiles
(
  const char   *name,       /*!<[in]      name of matching profile */
  int           n_coords,   /*!<[in]      number of point coordinates */
  cs_real_t     coords[][3] /*!<[in, out] point coordinates */
);

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function is used to compute user defined calculator formulae.
 *        The mathematical expression is defined in the GUI.
 */
/*----------------------------------------------------------------------------*/

void
cs_meg_post_calculator
(
  const char      *name,     /*!<[in] function name */
  const cs_lnum_t  n_elts,   /*!<[in] number of elements related to the zone */
  const cs_lnum_t *elt_ids,  /*!<[in] list of element ids related to the zone */
  const cs_real_t  xyz[][3], /*!<[in] list of coordinates related to the zone */
  cs_real_t       *retvals   /*!<[out] array of computed values */
);

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function is used to get the cs_cutcell_func_t pointer for a
 * given object.
 *
 * \return pointer to corresponding function.
 */
/*----------------------------------------------------------------------------*/

cs_cutcell_func_t *
cs_meg_ibm_func_by_name
(
  const char *object_name /*!<[in] name of the object */
);

/*----------------------------------------------------------------------------*/
/*!
 *
 * \brief This function is used to get the cs_ibm_volume_func_t pointer for a
 *        given object
 *
 * \return pointer to corresponding function.
 */
/*----------------------------------------------------------------------------*/

cs_ibm_volume_func_t *
cs_meg_ibm_volume_func_by_name
(
  const char *object_name, /*!<[in] Name of the immersed object */
  const char *gui_var_name /*!<[in] Name of the variable */
);

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 *
 * \brief This function is used to get the cs_ibm_fsi_func_t pointer for a
 *        given object
 *
 * \return pointer to corresponding function.
 */
/*----------------------------------------------------------------------------*/

cs_ibm_fsi_func_t *
cs_meg_ibm_fsi_func_by_name
(
  const char *object_name, /*!<[in] Name of the immersed object */
  const char *gui_var_name /*!<[in] Name of the variable */
);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MEG_PROTOTYPES_H__ */
